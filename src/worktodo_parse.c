#include "Mdata.h"

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static bool wp_fail(wp_cursor *c, const char *fmt, ...) {
	va_list ap;
	c->failed = true;
	va_start(ap, fmt);
	vsnprintf(c->err, sizeof(c->err), fmt, ap);
	va_end(ap);
	return false;
}

static void wp_skip_ws(wp_cursor *c) {
	while (c->s < c->limit && isspace((unsigned char)*c->s)) {
		c->s++;
	}
}

static bool wp_ieq_n(const char *s, const char *lit, size_t n) {
	size_t i;
	for (i = 0; i < n; i++) {
		if (tolower((unsigned char)s[i]) != tolower((unsigned char)lit[i])) {
			return false;
		}
	}
	return true;
}

void wp_cursor_init(wp_cursor *c, char *line) {
	memset(c, 0, sizeof(*c));
	c->s = line;
	c->limit = strchr(line, '\n');
	if (c->limit == NULL) {
		c->limit = line + strlen(line);
	}
	c->line_start = line;
	c->line_end = c->limit;
}

const char *wp_cursor_error(const wp_cursor *c) {
	return c->err;
}

bool wp_at_end(wp_cursor *c) {
	wp_cursor t = *c;
	wp_skip_ws(&t);
	return t.s >= t.limit;
}

bool wp_at_factors(wp_cursor *c) {
	wp_cursor t = *c;
	wp_skip_ws(&t);
	if (t.s >= t.limit || *t.s != ',') {
		return false;
	}
	t.s++;
	wp_skip_ws(&t);
	return t.s < t.limit && *t.s == '"';
}

bool wp_expect_comma(wp_cursor *c) {
	wp_skip_ws(c);
	if (c->s >= c->limit || *c->s != ',') {
		return wp_fail(c, "expected ','");
	}
	c->s++;
	return true;
}

bool wp_expect_uint(wp_cursor *c, uint64 *out) {
	char *end;
	wp_skip_ws(c);
	if (c->s >= c->limit || !isdigit((unsigned char)*c->s)) {
		return wp_fail(c, "expected an unsigned integer");
	}
	errno = 0;
	*out = strtoull(c->s, &end, 10);
	if (end == c->s || errno == ERANGE) {
		return wp_fail(c, "integer field out of range");
	}
	c->s = end;
	return true;
}

bool wp_expect_int(wp_cursor *c, int64 *out) {
	char *end;
	wp_skip_ws(c);
	if (c->s >= c->limit ||
	    !(isdigit((unsigned char)*c->s) ||
	      ((*c->s == '+' || *c->s == '-') && c->s + 1 < c->limit && isdigit((unsigned char)c->s[1])))) {
		return wp_fail(c, "expected an integer");
	}
	errno = 0;
	*out = strtoll(c->s, &end, 10);
	if (end == c->s || errno == ERANGE) {
		return wp_fail(c, "integer field out of range");
	}
	c->s = end;
	return true;
}

bool wp_expect_double(wp_cursor *c, double *out) {
	char *end;
	wp_skip_ws(c);
	if (c->s >= c->limit || !(isdigit((unsigned char)*c->s) || *c->s == '.')) {
		return wp_fail(c, "expected a number");
	}
	errno = 0;
	*out = strtod(c->s, &end);
	if (end == c->s || errno == ERANGE) {
		return wp_fail(c, "number field out of range");
	}
	c->s = end;
	return true;
}

bool wp_expect_quoted(wp_cursor *c, char **start, size_t *len) {
	char *q;
	wp_skip_ws(c);
	if (c->s >= c->limit || *c->s != '"') {
		return wp_fail(c, "expected a quoted string");
	}
	q = memchr(c->s + 1, '"', (size_t)(c->limit - c->s - 1));
	if (q == NULL) {
		return wp_fail(c, "unterminated quoted string");
	}
	if (start != NULL) {
		*start = c->s;
	}
	if (len != NULL) {
		*len = (size_t)(q - c->s) + 1;
	}
	c->s = q + 1;
	return true;
}

bool wp_try_double(wp_cursor *c, double *out) {
	wp_cursor save = *c;
	if (!wp_expect_comma(c)) {
		*c = save;
		return false;
	}
	wp_skip_ws(c);
	if (c->s >= c->limit || !(isdigit((unsigned char)*c->s) || *c->s == '.')) {
		*c = save;
		return false;
	}
	if (!wp_expect_double(c, out)) {
		*c = save;
		return false;
	}
	return true;
}

bool wp_try_uint(wp_cursor *c, uint64 *out) {
	wp_cursor save = *c;
	if (!wp_expect_comma(c)) {
		*c = save;
		return false;
	}
	wp_skip_ws(c);
	if (c->s >= c->limit || !isdigit((unsigned char)*c->s)) {
		*c = save;
		return false;
	}
	if (!wp_expect_uint(c, out)) {
		*c = save;
		return false;
	}
	return true;
}

/* k,b,n,c of the modulus m = k*b^n + c */
static bool wp_parse_kbnc(wp_cursor *c, wp_assignment *a) {
	return wp_expect_uint(c, &a->k) && wp_expect_comma(c) &&
	       wp_expect_uint(c, &a->b) && wp_expect_comma(c) &&
	       wp_expect_uint(c, &a->n) && wp_expect_comma(c) &&
	       wp_expect_int(c, &a->c);
}

/* ", \"f1,f2,...\"" -- a required trailing known-factors list */
static bool wp_take_factors(wp_cursor *c, wp_assignment *a) {
	if (!wp_expect_comma(c)) {
		return false;
	}
	wp_skip_ws(c);
	if (c->s >= c->limit || *c->s != '"') {
		return wp_fail(c, "expected a quoted known-factors list");
	}
	return wp_expect_quoted(c, &a->factors_start, &a->factors_len);
}

static bool wp_parse_optional_aid(wp_cursor *c, wp_assignment *a) {
	size_t i;
	wp_skip_ws(c);
	for (i = 0; i < WP_AID_MAX; i++) {
		if (c->s + i >= c->limit || !isxdigit((unsigned char)c->s[i])) {
			break;
		}
	}
	if (i == WP_AID_MAX && c->s + WP_AID_MAX < c->limit && c->s[WP_AID_MAX] == ',') {
		memcpy(a->aid, c->s, WP_AID_MAX);
		a->aid[WP_AID_MAX] = 0;
		a->has_aid = true;
		c->s += WP_AID_MAX + 1;
		return true;
	}
	if (c->s + 3 < c->limit && wp_ieq_n(c->s, "n/a", 3) && c->s[3] == ',') {
		c->s += 4;
		return true;
	}
	return true;
}

static bool wp_parse_keyword(wp_cursor *c, char *kw, size_t kwlen) {
	size_t n = 0;
	wp_skip_ws(c);
	while (c->s + n < c->limit && (isalnum((unsigned char)c->s[n]) || c->s[n] == '_')) {
		n++;
	}
	if (n == 0 || n + 1 > kwlen) {
		return wp_fail(c, "expected an assignment keyword");
	}
	memcpy(kw, c->s, n);
	kw[n] = 0;
	c->s += n;
	return true;
}

bool wp_parse_assignment(wp_cursor *c, wp_assignment *a) {
	char kw[32];
	uint64 v;

	memset(a, 0, sizeof(*a));

	/* Optional "Fermat," or "Mersenne," prefix (default Mersenne), then the keyword=,
	   then an optional 32-hexit PrimeNet assignment ID or "n/a" and a comma. */
	wp_skip_ws(c);
	if (c->s + 9 <= c->limit && wp_ieq_n(c->s, "Mersenne,", 9)) {
		a->mod_type = WP_MOD_MERSENNE;
		a->has_mod_prefix = true;
		c->s += 9;
	} else if (c->s + 7 <= c->limit && wp_ieq_n(c->s, "Fermat,", 7)) {
		a->mod_type = WP_MOD_FERMAT;
		a->has_mod_prefix = true;
		c->s += 7;
	} else {
		a->mod_type = WP_MOD_MERSENNE;
	}

	if (!wp_parse_keyword(c, kw, sizeof(kw))) {
		return false;
	}
	if (c->s >= c->limit || *c->s != '=') {
		return wp_fail(c, "expected '=' after '%s'", kw);
	}
	c->s++;

	if (!wp_parse_optional_aid(c, a)) {
		return false;
	}

	if (!strcasecmp(kw, "Test") || !strcasecmp(kw, "DoubleCheck")) {
		/* Test=p[,depth[,pminus1ed]]; DoubleCheck= has the same shape. */
		a->work = !strcasecmp(kw, "Test") ? WP_WORK_TEST : WP_WORK_DOUBLECHECK;
		if (!wp_expect_uint(c, &a->n)) return false;
		if (wp_try_double(c, &a->sieve_depth)) a->has_sieve_depth = true;
		/* Capture the pm1_done span so a P-1 split can rewrite just that field. */
		{
			wp_cursor probe = *c;
			char *field_start;
			wp_skip_ws(&probe);
			if (probe.s < probe.limit && *probe.s == ',') {
				probe.s++;
				wp_skip_ws(&probe);
				field_start = probe.s;
				if (wp_expect_uint(&probe, &v)) {
					a->pminus1ed = (int)v;
					a->has_pminus1ed = true;
					a->pminus1ed_pos = field_start;
					a->pminus1ed_end = probe.s;
					*c = probe;
				}
			}
		}
	} else if (!strcasecmp(kw, "AdvancedTest")) {
		/* AdvancedTest=p */
		a->work = WP_WORK_ADVANCEDTEST;
		if (!wp_expect_uint(c, &a->n)) return false;
	} else if (!strcasecmp(kw, "Factor")) {
		/* Factor=p[,depth[,factor_to]] */
		a->work = WP_WORK_FACTOR;
		if (!wp_expect_uint(c, &a->n)) return false;
		if (wp_try_double(c, &a->sieve_depth)) a->has_sieve_depth = true;
		if (wp_try_double(c, &a->factor_to)) a->has_factor_to = true;
	} else if (!strcasecmp(kw, "Pfactor")) {
		/* Pfactor=k,b,n,c,depth,tests_saved[,"factors"] */
		a->work = WP_WORK_PFACTOR;
		if (!wp_parse_kbnc(c, a)) return false;
		if (!wp_expect_comma(c) || !wp_expect_double(c, &a->sieve_depth)) return false;
		a->has_sieve_depth = true;
		if (!wp_expect_comma(c) || !wp_expect_double(c, &a->tests_saved)) return false;
		a->has_tests_saved = true;
		if (wp_at_factors(c) && !wp_take_factors(c, a)) return false;
	} else if (!strcasecmp(kw, "Pminus1")) {
		/* Pminus1=k,b,n,c,B1,B2[,depth[,B2_start]][,"factors"] */
		a->work = WP_WORK_PMINUS1;
		if (!wp_parse_kbnc(c, a)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &a->b1)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &a->b2)) return false;
		if (wp_try_double(c, &a->sieve_depth)) a->has_sieve_depth = true;
		if (wp_try_uint(c, &v)) { a->b2_start = v; a->has_b2_start = (v > a->b1); }
		if (wp_at_factors(c) && !wp_take_factors(c, a)) return false;
	} else if (!strcasecmp(kw, "Pplus1")) {
		/* Pplus1=k,b,n,c,B1,B2,nth_run[,depth][,"factors"] */
		a->work = WP_WORK_PPLUS1;
		if (!wp_parse_kbnc(c, a)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &a->b1)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &a->b2)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &v)) return false;	/* nth_run */
		if (wp_try_double(c, &a->sieve_depth)) a->has_sieve_depth = true;
		if (wp_at_factors(c) && !wp_take_factors(c, a)) return false;
	} else if (!strcasecmp(kw, "PRP") || !strcasecmp(kw, "PRPDC")) {
		/* PRP=k,b,n,c[,depth,tests_saved[,base,residue_type]][,"factors"]; PrimeNet's
		   PRPDC= double-check assignments use the same shape. */
		a->work = !strcasecmp(kw, "PRP") ? WP_WORK_PRP : WP_WORK_PRPDC;
		if (!wp_parse_kbnc(c, a)) return false;
		if (wp_at_factors(c)) {
			if (!wp_take_factors(c, a)) return false;
		} else if (!wp_at_end(c)) {
			if (!wp_expect_comma(c) || !wp_expect_double(c, &a->sieve_depth)) return false;
			a->has_sieve_depth = true;
			if (!wp_expect_comma(c)) return false;
			wp_skip_ws(c);
			a->tests_saved_pos = c->s;
			if (!wp_expect_double(c, &a->tests_saved)) return false;
			a->tests_saved_end = c->s;
			a->has_tests_saved = true;
			if (wp_at_factors(c)) {
				if (!wp_take_factors(c, a)) return false;
			} else if (!wp_at_end(c)) {
				int64 rt;
				if (!wp_expect_comma(c) || !wp_expect_uint(c, &v)) return false;
				a->prp_base = (uint32)v;
				if (!wp_expect_comma(c) || !wp_expect_int(c, &rt)) return false;
				a->residue_type = (uint32)rt;
				a->has_prp_pair = true;
				if (wp_at_factors(c) && !wp_take_factors(c, a)) return false;
			}
		}
	} else if (!strcasecmp(kw, "ECM") || !strcasecmp(kw, "ECM2")) {
		/* ECM=k,b,n,c,B1[,additional numbers][,"factors"] */
		a->work = WP_WORK_ECM;
		if (!wp_parse_kbnc(c, a)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &a->b1)) return false;
		while (!wp_at_end(c) && !wp_at_factors(c)) {
			if (!wp_expect_comma(c) || !wp_expect_uint(c, &v)) return false;
		}
		if (wp_at_factors(c) && !wp_take_factors(c, a)) return false;
	} else if (!strcasecmp(kw, "Cert")) {
		/* Cert=k,b,n,c,num_squarings */
		a->work = WP_WORK_CERT;
		if (!wp_parse_kbnc(c, a)) return false;
		if (!wp_expect_comma(c) || !wp_expect_uint(c, &v)) return false;	/* num_squarings */
	} else {
		return wp_fail(c, "unrecognized worktodo keyword '%s'", kw);
	}

	if (!wp_at_end(c)) {
		return wp_fail(c, "trailing junk after the assignment fields");
	}
	return true;
}

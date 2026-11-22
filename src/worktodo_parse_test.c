#include "Mdata.h"

#include <stdio.h>
#include <string.h>

static int failures;
static wp_assignment a;
static char err[WP_ERRLEN];

static int parse_line(const char *line) {
	/* Static: wp_parse_assignment() stores pointers into it, read by the assertions after return. */
	static char buf[1024];
	wp_cursor c;

	strncpy(buf, line, sizeof(buf) - 1);
	buf[sizeof(buf) - 1] = 0;
	wp_cursor_init(&c, buf);
	if (!wp_parse_assignment(&c, &a)) {
		snprintf(err, sizeof(err), "%s", wp_cursor_error(&c));
		return 0;
	}
	return 1;
}

#define T(cond) do { if (!(cond)) { printf("FAIL [%s]: %s\n", name, #cond); failures++; } } while (0)
#define CASE(str) { const char *name = (str); \
	if (!parse_line(name)) { printf("FAIL [%s]: parse error: %s\n", name, err); failures++; } else {

#define ENDCASE } }

/* A line that must be rejected: the parse has to fail. */
#define FAILCASE(str) { const char *name = (str); \
	if (parse_line(name)) { printf("FAIL [%s]: expected a parse error\n", name); failures++; } }

int main(void) {
	CASE("Test=12345")
		T(a.work == WP_WORK_TEST);
		T(a.mod_type == WP_MOD_MERSENNE);
		T(a.n == 12345);
		T(!a.has_sieve_depth && !a.has_pminus1ed);
	ENDCASE;

	CASE("Test=12345,75,0")
		T(a.work == WP_WORK_TEST);
		T(a.n == 12345);
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_pminus1ed && a.pminus1ed == 0);
	ENDCASE;

	CASE("Test=12345 ,75,0")
		T(a.n == 12345);
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_pminus1ed && a.pminus1ed == 0);
	ENDCASE;

	CASE("DoubleCheck=12345,70,1")
		T(a.work == WP_WORK_DOUBLECHECK);
		T(a.n == 12345);
		T(a.has_sieve_depth && a.sieve_depth == 70.0);
		T(a.has_pminus1ed && a.pminus1ed == 1);
	ENDCASE;

	CASE("Factor=12345 ,70,75")
		T(a.work == WP_WORK_FACTOR);
		T(a.n == 12345);
		T(a.has_sieve_depth && a.sieve_depth == 70.0);
		T(a.has_factor_to && a.factor_to == 75.0);
	ENDCASE;

	CASE("Factor=12345,70,75")
		T(a.work == WP_WORK_FACTOR);
		T(a.n == 12345);
		T(a.has_sieve_depth && a.sieve_depth == 70.0);
		T(a.has_factor_to && a.factor_to == 75.0);
	ENDCASE;

	CASE("PRP=1,2,2063,-1")
		T(a.work == WP_WORK_PRP);
		T(a.k == 1 && a.b == 2 && a.n == 2063 && a.c == -1);
		T(!a.has_sieve_depth && !a.has_tests_saved && !a.has_prp_pair);
		T(a.factors_start == NULL);
	ENDCASE;

	CASE("PRP=C42540C352E54E906108D48FA5D89488,1,2,80340397,-1,75,0,3,1")
		T(a.work == WP_WORK_PRP);
		T(a.has_aid && strcmp(a.aid, "C42540C352E54E906108D48FA5D89488") == 0);
		T(a.n == 80340397 && a.c == -1);
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 0.0);
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 1);
		T(a.factors_start == NULL);
	ENDCASE;

	CASE("PRP=1,2,3215747,-1,99,0,3,5,\"4457025343,185822885311153245017\"")
		T(a.n == 3215747);
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 5);
		T(a.factors_start != NULL);
		T(a.factors_len == strlen("\"4457025343,185822885311153245017\""));
		T(strncmp(a.factors_start, "\"4457025343,185822885311153245017\"", a.factors_len) == 0);
	ENDCASE;

	CASE("PRP=n/a,1,2,33554432,+1,99,0,3,5,\"2170072644496392193\"")
		T(!a.has_aid);
		T(a.n == 33554432 && a.c == 1);
		T(a.has_prp_pair && a.residue_type == 5);
		T(a.factors_len == strlen("\"2170072644496392193\""));
	ENDCASE;

	CASE("PRP=1,2,2063,-1, 75,0,3,1")
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 0.0);
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 1);
	ENDCASE;

	CASE("PRP=1,2,2063,-1,\"4127\"")
		T(!a.has_sieve_depth && !a.has_tests_saved && !a.has_prp_pair);
		T(a.factors_start != NULL && a.factors_len == 6);
		T(strncmp(a.factors_start, "\"4127\"", 6) == 0);
	ENDCASE;

	CASE("PRP=1,2,2063,-1,99,0,3,5,\"4127\"   \r\n")
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 5);
		T(a.factors_start != NULL && a.factors_len == 6);
	ENDCASE;

	/* The PR #55 acceptance matrix from MarkRose, with the AID dropped and p = M120371.
	   Grammar only: the [0,2] tests_saved range is enforced in the dispatch code, not here. */
	CASE("PRP=1,2,120371,-1")
		T(a.work == WP_WORK_PRP && a.n == 120371);
		T(!a.has_sieve_depth && !a.has_tests_saved && !a.has_prp_pair);
		T(a.factors_start == NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,75,0")
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 0.0);
		T(a.factors_start == NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,75,0,3,1")
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 1);
		T(a.factors_start == NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1, 75,0,3,1")
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 0.0);
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 1);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,\"240743\"")
		T(!a.has_sieve_depth && !a.has_tests_saved && !a.has_prp_pair);
		T(a.factors_start != NULL && strstr(a.factors_start, "240743") != NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,\"240743,152389687\"\r\n")
		T(!a.has_sieve_depth && !a.has_tests_saved);
		T(a.factors_start != NULL && strstr(a.factors_start, "152389687") != NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,75,0,\"240743\"")
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 0.0);
		T(a.factors_start != NULL && strstr(a.factors_start, "240743") != NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,75,0,3,5,\"240743\"")
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 5);
		T(a.factors_start != NULL && strstr(a.factors_start, "240743") != NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,75,1")
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 1.0);
		T(!a.has_prp_pair && a.factors_start == NULL);
	ENDCASE;

	CASE("PRP=1,2,120371,-1,75,1,3,5,\"240743\"")
		T(a.has_tests_saved && a.tests_saved == 1.0);
		T(a.has_prp_pair && a.prp_base == 3 && a.residue_type == 5);
		T(a.factors_start != NULL && strstr(a.factors_start, "240743") != NULL);
	ENDCASE;

	/* Out of range for the dispatch code, but syntactically fine: the parser reports the value. */
	CASE("PRP=1,2,120371,-1,75,5")
		T(a.has_tests_saved && a.tests_saved == 5.0);
	ENDCASE;

	CASE("Pminus1=1,2,2063,-1,1000,10000")
		T(a.work == WP_WORK_PMINUS1);
		T(a.n == 2063 && a.c == -1);
		T(a.b1 == 1000 && a.b2 == 10000);
		T(!a.has_sieve_depth && !a.has_b2_start);
	ENDCASE;

	CASE("Pminus1=1,2,2063,-1,1000,10000,75,2000,\"4127\"")
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_b2_start && a.b2_start == 2000);
		T(a.factors_start != NULL && a.factors_len == 6);
	ENDCASE;

	CASE("Pminus1=1,2,2063,-1,1000,10000,\"4127\"")
		T(a.b1 == 1000 && a.b2 == 10000);
		T(!a.has_sieve_depth && !a.has_b2_start);
		T(a.factors_start != NULL && a.factors_len == 6);
	ENDCASE;

	CASE("Pfactor=1,2,2063,-1,75,1.3")
		T(a.work == WP_WORK_PFACTOR);
		T(a.has_sieve_depth && a.sieve_depth == 75.0);
		T(a.has_tests_saved && a.tests_saved == 1.3);
	ENDCASE;

	CASE("Pfactor=1,2,2063,-1,75,1.3,\"4127\"")
		T(a.has_tests_saved && a.tests_saved == 1.3);
		T(a.factors_start != NULL && a.factors_len == 6);
	ENDCASE;

	CASE("Fermat,Test=25")
		T(a.work == WP_WORK_TEST && a.mod_type == WP_MOD_FERMAT && a.has_mod_prefix);
		T(a.n == 25);
	ENDCASE;

	CASE("Mersenne,Test=12345")
		T(a.work == WP_WORK_TEST && a.mod_type == WP_MOD_MERSENNE && a.has_mod_prefix);
		T(a.n == 12345);
	ENDCASE;

	CASE("PRPDC=1,2,2063,-1")
		T(a.work == WP_WORK_PRPDC);
		T(a.n == 2063);
	ENDCASE;

	CASE("ECM=1,2,2063,-1,11000,0,100,12345,\"4127\"")
		T(a.work == WP_WORK_ECM);
		T(a.b1 == 11000);
		T(a.factors_start != NULL);
	ENDCASE;

	CASE("Cert=1,2,2063,-1,12345")
		T(a.work == WP_WORK_CERT);
		T(a.n == 2063);
	ENDCASE;

	FAILCASE("Foo=1");
	FAILCASE("PRP=1,2,2063,-1,99");
	FAILCASE("PRP=1,2,2063,-1,99,0,3");
	FAILCASE("PRP=1,2,2063,-1,99,0,3,1,5");

	printf("%d failures\n", failures);
	return failures != 0;
}

////////////////////////////////////////////
//
// secure5.c #include module
//
// PrimeNet v5 client API security module
// Kurowski 09/28/2005
// Woltman 06/13/2007 - use md5 code already in prime95, eliminated globals
//
// 32-bit UNIQUE TRUST CLIENT BUILD CONSTANT
// Randomly choose an available 32-bit number for this constant,
// which is set in the #define below before this source module
// is given to a trusted client application builder.
// DO NOT SHARE CONSTANTS among trusted builders!
#define _V5_UNIQUE_TRUSTED_CLIENT_CONSTANT_	17737	// GIMPS WOLTMAN
////////////////////////////////////////////

#include <time.h>

#define _V5_SECURITY_MODULE_PRESENT_

void make_v5_client_key (char *p_v5key, char computer_guid[33])
{
	unsigned char k[16];
	int i;

	md5_raw_output (k, computer_guid);
	for (i = 0; i < 16; i++)
		k[i] ^= k[(k[i] ^ (_V5_UNIQUE_TRUSTED_CLIENT_CONSTANT_ & 0xFF)) % 16]
			 ^ (_V5_UNIQUE_TRUSTED_CLIENT_CONSTANT_ / 256);
	md5_raw_input (p_v5key, k, 16);
	strupper (p_v5key);
}

// converts URL to secured v5_URL

void secure_v5_url (char *URL, char *p_v5key)
{
	char hash[33];
	unsigned int len;
	static unsigned int __v5salt_ = 0;

	len = (unsigned int) strlen (URL);
	if (len > 4020) {
		URL[0] = 0;	// return empty URL if URL is too long
		return;
	}

	if (!__v5salt_) srand ((unsigned int) time(0));
	__v5salt_ = rand () & 0xFFFF;

	sprintf (URL + len, "&ss=%u&%s", __v5salt_, p_v5key);

	md5 (hash, URL);
	strupper (hash);
	sprintf (URL + len, "&ss=%u&sh=%s", __v5salt_, hash);
}

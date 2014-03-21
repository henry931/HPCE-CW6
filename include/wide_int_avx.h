#ifndef wide_int_avx_h
#define wide_int_avx_h

#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

// AVX Header
#include "immintrin.h"

int wide_compare_avx(__m256i *a, __m256i *b)
{
	if (a == b)
		return 0;

	uint32_t* const a_ptr = (uint32_t*)a;
	uint32_t* const b_ptr = (uint32_t*)b;

	if (a_ptr[7]>b_ptr[7]) return +1;
	if (a_ptr[7]<b_ptr[7]) return -1;
	if (a_ptr[6]>b_ptr[6]) return +1;
	if (a_ptr[6]<b_ptr[6]) return -1;
	if (a_ptr[5]>b_ptr[5]) return +1;
	if (a_ptr[5]<b_ptr[5]) return -1;
	if (a_ptr[4]>b_ptr[4]) return +1;
	if (a_ptr[4]<b_ptr[4]) return -1;
	if (a_ptr[3]>b_ptr[3]) return +1;
	if (a_ptr[3]<b_ptr[3]) return -1;
	if (a_ptr[2]>b_ptr[2]) return +1;
	if (a_ptr[2]<b_ptr[2]) return -1;
	if (a_ptr[1]>b_ptr[1]) return +1;
	if (a_ptr[1]<b_ptr[1]) return -1;
	if (a_ptr[0]>b_ptr[0]) return +1;
	if (a_ptr[0]<b_ptr[0]) return -1;

	return 0;
}

void wide_xor_avx(__m256i *result, __m256i *a, __m256i *b)
{
	*result = *(__m256i*)(&_mm256_xor_ps(*(__m256*)a, *(__m256*)b));
}

void wide_zero_avx(__m256i *result)
{
	*result = _mm256_setzero_si256();
}

void wide_ones_avx(__m256i *result)
{
	*result = _mm256_set1_epi32(0xFFFFFFFFul);
}

void wide_copy_avx(__m256i *result, __m256i *a)
{
	*result = *a;
}

void wide_mul_avx(__m256i *result, __m256i *a, __m256i *b)
{
	// So the idea is to do long multiplication...

	// Pack 32-bit limbs
	__m128i* const ainput = (__m128i*) a;
	__m128i* const binput = (__m128i*) b;
	const __m128i aupper = _mm_cvtepu32_epi64(_mm_srli_si128(*ainput, 8));
	const __m128i alower = _mm_cvtepu32_epi64(*ainput);
	const __m128i bupper_shift0 = _mm_cvtepu32_epi64(_mm_srli_si128(*binput, 8));	// also blower_shift2
	const __m128i blower_shift0 = _mm_cvtepu32_epi64(*binput);						// also bupper_shift2
	const __m128i bupper_shift1 = _mm_cvtepu32_epi64(_mm_shuffle_epi32(*binput, _MM_SHUFFLE(0, 0, 0, 3))); // also blower_shift3
	const __m128i blower_shift1 = _mm_cvtepu32_epi64(_mm_shuffle_epi32(*binput, _MM_SHUFFLE(0, 0, 2, 1))); // also bupper_shift3

	// Multiply limbs
	const __m128i abupper_shift0 = _mm_mul_epu32(aupper, bupper_shift0);
	const __m128i ablower_shift0 = _mm_mul_epu32(alower, blower_shift0);
	const __m128i abupper_shift1 = _mm_mul_epu32(aupper, bupper_shift1);
	const __m128i ablower_shift1 = _mm_mul_epu32(alower, blower_shift1);
	const __m128i abupper_shift2 = _mm_mul_epu32(aupper, blower_shift0);
	const __m128i ablower_shift2 = _mm_mul_epu32(alower, bupper_shift0);
	const __m128i abupper_shift3 = _mm_mul_epu32(aupper, blower_shift1);
	const __m128i ablower_shift3 = _mm_mul_epu32(alower, bupper_shift1);

	// Convert to 64-bit for addition overflow
	const __m128i r0c0 = _mm_cvtepu32_epi64(ablower_shift0);
	const __m128i r0c1 = _mm_cvtepu32_epi64(ablower_shift1);
	const __m128i r0c2 = _mm_cvtepu32_epi64(ablower_shift2);
	const __m128i r2c2 = _mm_cvtepu32_epi64(abupper_shift2);
	const __m128i r0c3 = _mm_cvtepu32_epi64(ablower_shift3);
	const __m128i r2c3 = _mm_cvtepu32_epi64(abupper_shift3);
	const __m128i r2c4 = _mm_cvtepu32_epi64(abupper_shift0);
	const __m128i r2c5 = _mm_cvtepu32_epi64(abupper_shift1);

	// Need to shift these guys 8-bytes right first for the rest
	const __m128i r1c1 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift3, 8));
	const __m128i r1c2 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift0, 8));
	const __m128i r1c3 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift1, 8));
	const __m128i r3c3 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift1, 8));
	const __m128i r1c4 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift2, 8));
	const __m128i r3c4 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift2, 8));
	const __m128i r3c5 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift3, 8));
	const __m128i r3c6 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift0, 8));

	// OK now add
	__m128i result0 = r0c0;
	__m128i result1 = _mm_add_epi64(r0c1, r1c1);
	__m128i result2 = _mm_add_epi64(_mm_add_epi64(r0c2, r1c2), r2c2);
	__m128i result3 = _mm_add_epi64(_mm_add_epi64(r0c3, r1c3), _mm_add_epi64(r2c3, r3c3));
	__m128i result4 = _mm_add_epi64(_mm_add_epi64(r1c4, r2c4), r3c4);
	__m128i result5 = _mm_add_epi64(r2c5, r3c5);
	__m128i result6 = r3c6;

	// Pointers for access
	uint32_t* const res0 = (uint32_t*)&result0;
	uint32_t* const res1 = (uint32_t*)&result1;
	uint32_t* const res2 = (uint32_t*)&result2;
	uint32_t* const res3 = (uint32_t*)&result3;
	uint32_t* const res4 = (uint32_t*)&result4;
	uint32_t* const res5 = (uint32_t*)&result5;
	uint32_t* const res6 = (uint32_t*)&result6;

	uint64_t* const res0_64 = (uint64_t*)&result0;
	uint64_t* const res1_64 = (uint64_t*)&result1;
	uint64_t* const res2_64 = (uint64_t*)&result2;
	uint64_t* const res3_64 = (uint64_t*)&result3;
	uint64_t* const res4_64 = (uint64_t*)&result4;
	uint64_t* const res5_64 = (uint64_t*)&result5;
	uint64_t* const res6_64 = (uint64_t*)&result6;

	// Push overflows into subsequent buffers
	res1_64[0] += res0_64[1] + res0[1];
	res2_64[0] += res1_64[1] + res1[1];
	res3_64[0] += res2_64[1] + res2[1];
	res4_64[0] += res3_64[1] + res3[1];
	res5_64[0] += res4_64[1] + res4[1];
	res6_64[0] += res5_64[1] + res5[1];

	// Output lowest 32-bit in each buffer
	*result = _mm256_set_epi32(res6[2] + res6[1], res6[0], res5[0], res4[0], res3[0], res2[0], res1[0], res0[0]);
}

double wide_as_double_avx(__m256i *a)
{
	uint32_t* const a_ptr = (uint32_t*)a;

	double acc = (double)a_ptr[0];
	acc += ldexp((double)a_ptr[1], 1 * 32);
	acc += ldexp((double)a_ptr[2], 2 * 32);
	acc += ldexp((double)a_ptr[3], 3 * 32);
	acc += ldexp((double)a_ptr[4], 4 * 32);
	acc += ldexp((double)a_ptr[5], 5 * 32);
	acc += ldexp((double)a_ptr[6], 6 * 32);
	acc += ldexp((double)a_ptr[7], 7 * 32);

	return acc;
}

// 2 * 128-bit input (crushed 256) -> 128-bit output + carry
uint32_t wide_add_avx_4limb(__m256i *result, __m256i *a, __m256i *b)
{
	const __m128i a_upper = _mm_cvtepu32_epi64(_mm_srli_si128(*(__m128i*)a, 8));
	const __m128i a_lower = _mm_cvtepu32_epi64(*(__m128i*)a);
	const __m128i b_upper = _mm_cvtepu32_epi64(_mm_srli_si128(*(__m128i*)b, 8));
	const __m128i b_lower = _mm_cvtepu32_epi64(*(__m128i*)b);
	__m128i result_upper = _mm_add_epi64(a_upper, b_upper);
	__m128i result_lower = _mm_add_epi64(a_lower, b_lower);
	uint32_t* const result_upper_ptr = (uint32_t*)&result_upper;
	uint32_t* const result_lower_ptr = (uint32_t*)&result_lower;
	uint64_t* const result_upper_ptr_64 = (uint64_t*)&result_upper;
	uint64_t* const result_lower_ptr_64 = (uint64_t*)&result_lower;

	result_lower_ptr_64[1] += result_lower_ptr[1];
	result_upper_ptr_64[0] += result_lower_ptr[3];
	result_upper_ptr_64[1] += result_upper_ptr[1];

	*result = _mm256_set_epi32(0, 0, 0, 0, result_upper_ptr[2], result_upper_ptr[0], result_lower_ptr[2], result_lower_ptr[0]);
	return result_upper_ptr[3];
}

// Just input whole vectors - none of this +4 nonsense
void PoolHashStep_wide_add(__m256i *x, __m256i *tmp)
{
	__m128i* const xsplit_ptr = (__m128i*)x;
	__m128i* const tmpsplit_ptr = (__m128i*)tmp;
	const __m128i a_upper = _mm_cvtepu32_epi64(_mm_srli_si128(tmpsplit_ptr[0], 8));
	const __m128i a_lower = _mm_cvtepu32_epi64(tmpsplit_ptr[0]);
	const __m128i b_upper = _mm_cvtepu32_epi64(_mm_srli_si128(xsplit_ptr[1], 8));
	const __m128i b_lower = _mm_cvtepu32_epi64(xsplit_ptr[1]);

	__m128i result_upper = _mm_add_epi64(a_upper, b_upper);
	__m128i result_lower = _mm_add_epi64(a_lower, b_lower);
	uint32_t* const result_upper_ptr = (uint32_t*)&result_upper;
	uint32_t* const result_lower_ptr = (uint32_t*)&result_lower;
	uint64_t* const result_upper_ptr_64 = (uint64_t*)&result_upper;
	uint64_t* const result_lower_ptr_64 = (uint64_t*)&result_lower;

	result_lower_ptr_64[1] += result_lower_ptr[1];
	result_upper_ptr_64[0] += result_lower_ptr[3];
	result_upper_ptr_64[1] += result_upper_ptr[1];

	uint32_t* const tmp_ptr = (uint32_t*)tmp;

	uint64_t acc = tmp_ptr[4] + result_upper_ptr[3];
	const uint32_t result4 = (uint32_t)acc;
	acc = (acc >> 32) + tmp_ptr[5];
	const uint32_t result5 = (uint32_t)acc;
	acc = (acc >> 32) + tmp_ptr[6];
	const uint32_t result6 = (uint32_t)acc;
	acc = (acc >> 32) + tmp_ptr[7];
	const uint32_t result7 = (uint32_t)acc;

	*x = _mm256_set_epi32(result7, result6, result5, result4, result_upper_ptr[2], result_upper_ptr[0], result_lower_ptr[2], result_lower_ptr[0]);
}

// 1 * 256-bit input, 64-bit chainHash, salt, Id, index
void PoolHash_wide_add(__m256i *InA, uint64_t *InB3, uint64_t *InB2, uint64_t *InB1, uint64_t *InB0)
{
	*InA = _mm256_set_epi32((uint32_t)*InB3, (uint32_t)*InB3, (uint32_t)*InB2, (uint32_t)*InB2, (uint32_t)*InB1, (uint32_t)*InB1, (uint32_t)*InB0, (uint32_t)*InB0);
}

#endif

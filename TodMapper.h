#include <math.h>
#include "bsalign/filereader.h"

void u8iseq_to_u1iseq(int len1, u8i *seq1, int len2, u1i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u8i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void u4iseq_to_u1iseq(int len1, u4i *seq1, int len2, u1i *seq2, int len)
{
	int i, j, p;
	if(len1 < len2)
	{
		p = (u8i)pow(2, len1);
		for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
	}
	else
	{
		p = (u8i)pow(2, len2);
		u4i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
	}
}

void u2iseq_to_u1iseq(int len1, u2i *seq1, int len2, u1i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u2i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void u1iseq_to_u1iseq(int len1, u1i *seq1, int len2, u1i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u1i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void u1iseq_to_u2iseq(int len1, u1i *seq1, int len2, u2i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u1i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void u1iseq_to_u4iseq(int len1, u1i *seq1, int len2, u4i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u1i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void u1iseq_to_u8iseq(int len1, u1i *seq1, int len2, u8i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u1i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void intseq_to_u1iseq(int len1, int *seq1, int len2, u1i *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                int tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

void u1iseq_to_intseq(int len1, u1i *seq1, int len2, int *seq2, int len)
{
        int i, j, p;
        if(len1 < len2)
        {
                p = (u8i)pow(2, len1);
                for(i = 0; i < len; i++)
                {
                        seq2[i] = 0;
                        for(j = 0; j < len2 / len1; j++)
                        {
                                seq2[i] = seq2[i] + seq1[i * (len2 / len1) + j] * (u8i)pow(p, (len2 / len1) - j - 1);
                        }
                }
        }
        else
        {
                p = (u8i)pow(2, len2);
                u1i tmp1, tmp2;
                for(i = 0; i < len; i++)
                {
                        tmp1 = seq1[i];
                        for(j = 0; j < len1 / len2; j++)
                        {
                                tmp2 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp2 = tmp1 - tmp2;
                                seq2[i * (len1 / len2) + j] = tmp2 / (u8i)pow(p, (len1 / len2) - j - 1);
                                tmp1 = tmp1 % (u8i)pow(p, (len1 / len2) - j - 1);
                        }
                }
        }
}

uint64_t hash64(uint64_t key, uint64_t mask)
{
        key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
        key = key ^ key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
}

struct MINIMIZER_T
{
	u8i s;
        u8i ip;
};

struct MINIMIZER_V
{
        u8i size;
	u8i cap;
	struct MINIMIZER_T *buffer;
};

void minimizer_push(struct MINIMIZER_V *v, struct MINIMIZER_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct MINIMIZER_T *)realloc(v->buffer, sizeof(struct MINIMIZER_T) * v->cap);
	}
	v->buffer[v->size++] = t;
}

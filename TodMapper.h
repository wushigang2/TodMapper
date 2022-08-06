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

struct QRYBIN_T
{
	u8i s;
};

struct QRYBIN_V
{
	u8i size;
	u8i cap;
	struct QRYBIN_T *buffer;
};

void qrybin_push(struct QRYBIN_V *v, struct QRYBIN_T t)
{
	u8i ui;
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct QRYBIN_T *)realloc(v->buffer, sizeof(struct QRYBIN_T) * v->cap);
	}
	for(ui = 0; ui < v->size; ui++)
	{
		if(v->buffer[ui].s == t.s)
		{
			v->buffer[ui] = v->buffer[v->size - 1];
			v->size--;
			break;
		}
	}
	v->buffer[v->size++] = t;
}

void getqrybin(BioSequence *seq, b4i ksize, b4i wsize, b4i bsize, b4i bi, struct MINIMIZER_V *minimizer_v, struct QRYBIN_V **qrybin_v)
{
	b4i pre_idx, min_idx, bj, bk, bl;
	u1i nucleobase[256];
	u8i key[4], mask[4], ui;
	struct MINIMIZER_T pre_elm[64], min_elm;
	struct QRYBIN_T qrybin_t;
	//
	for(bj = 0; bj < 128; bj++)
	{
		nucleobase[bj] = 4;
	}
	nucleobase[65] = nucleobase[97] = 0;
	nucleobase[67] = nucleobase[99] = 1;
	nucleobase[71] = nucleobase[103] = 2;
	nucleobase[84] = nucleobase[85] = nucleobase[116] = nucleobase[117] = 3;
	mask[0] = (1ULL << (ksize * 2)) - 1;
	mask[1] = (1ULL << 62) - 1;
	mask[2] = (ksize - 1) * 2;
	mask[3] = (1ULL << 32) - 1;
	//
	pre_idx = min_idx = 0;
	key[0] = key[1] = 0;
	min_elm.s = UINT64_MAX;
	bk = 0;
	for(bj = 0; bj < seq->seq->size; bj++)
	{
		nucleobase[128] = nucleobase[(int)seq->seq->string[bj]];
		if(nucleobase[128] < 4)
		{
			key[0] = ((key[0] << 2) & mask[0]) | nucleobase[128];
			key[1] = ((key[1] >> 2) & mask[1]) | ((nucleobase[128] ^ 3ULL) << mask[2]);
			bk++;
			if(bk >= ksize)
			{
				key[2] = hash64(key[0], mask[0]);
				key[3] = hash64(key[1], mask[0]);
				pre_elm[pre_idx].s = key[2] <= key[3] ? key[2] : key[3];
				pre_elm[pre_idx].ip = ((u8i)bi << 32) | (bj - ksize + 1);
				minimizer_push(minimizer_v, pre_elm[pre_idx]);
				if(bk < ksize + wsize - 1)
				{
					if(min_elm.s >= pre_elm[pre_idx].s)
					{
						min_elm = pre_elm[pre_idx];
						min_idx = pre_idx;
					}
					pre_idx = pre_idx == wsize - 1 ? 0 : pre_idx + 1;
				}
				else if(bk == ksize + wsize - 1)
				{
					if(min_elm.s >= pre_elm[pre_idx].s)
					{
						min_elm = pre_elm[pre_idx];
						min_idx = pre_idx;
					}
					pre_idx = pre_idx == wsize - 1 ? 0 : pre_idx + 1;
					for(bl = pre_idx; bl < wsize; bl++)
					{
						if(pre_elm[bl].s == min_elm.s)
						{
							qrybin_t.s = pre_elm[bl].s;
							ui = pre_elm[bl].ip & mask[3];
							ui = ui - (ui % bsize);
							ui = ui / bsize;
							qrybin_push(qrybin_v[ui], qrybin_t);
						}
					}
					for(bl = 0; bl < pre_idx; bl++)
					{
						if(pre_elm[bl].s == min_elm.s)
						{
							qrybin_t.s = pre_elm[bl].s;
							ui = pre_elm[bl].ip & mask[3];
							ui = ui - (ui % bsize);
							ui = ui / bsize;
							qrybin_push(qrybin_v[ui], qrybin_t);
						}
					}
				}
				else
				{
					if(min_elm.s >= pre_elm[pre_idx].s)
					{
						qrybin_t.s = pre_elm[pre_idx].s;
						ui = pre_elm[pre_idx].ip & mask[3];
						ui = ui - (ui % bsize);
						ui = ui / bsize;
						qrybin_push(qrybin_v[ui], qrybin_t);
						min_elm = pre_elm[pre_idx];
						min_idx = pre_idx;
						pre_idx = pre_idx == wsize - 1 ? 0 : pre_idx + 1;
					}
					else
					{
						if(min_idx == pre_idx)
						{
							min_elm.s = UINT64_MAX;
							pre_idx = pre_idx == wsize - 1 ? 0 : pre_idx + 1;
							for(bl = pre_idx; bl < wsize; bl++)
							{
								if(min_elm.s >= pre_elm[bl].s)
								{
									min_elm = pre_elm[bl];
									min_idx = bl;
								}
							}
							for(bl = 0; bl < pre_idx; bl++)
							{
								if(min_elm.s >= pre_elm[bl].s)
								{
									min_elm = pre_elm[bl];
									min_idx = bl;
								}
							}
							for(bl = pre_idx; bl < wsize; bl++)
							{
								if(pre_elm[bl].s == min_elm.s)
								{
									qrybin_t.s = pre_elm[bl].s;
									ui = pre_elm[bl].ip & mask[3];
									ui = ui - (ui % bsize);
									ui = ui / bsize;
									qrybin_push(qrybin_v[ui], qrybin_t);
								}
							}
							for(bl = 0; bl < pre_idx; bl++)
							{
								if(pre_elm[bl].s == min_elm.s)
								{
									qrybin_t.s = pre_elm[bl].s;
									ui = pre_elm[bl].ip & mask[3];
									ui = ui - (ui % bsize);
									ui = ui / bsize;
									qrybin_push(qrybin_v[ui], qrybin_t);
								}
							}
						}
						else
						{
							pre_idx = pre_idx == wsize - 1 ? 0 : pre_idx + 1;
						}
					}
				}
			}
			else
			{
			}
		}
		else
		{
			bk = 0;
			min_elm.s = UINT64_MAX;
		}
	}
}

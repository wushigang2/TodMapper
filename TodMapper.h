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
			break;
		}
	}
	if(ui == v->size)
	{
		v->buffer[v->size++] = t;
	}
}

struct ANCHOR_T
{
	b4i x;
	b4i y;
	float f;
	b4i p;
};

struct ANCHOR_V
{
	u8i size;
	u8i cap;
	struct ANCHOR_T *buffer;
};

void anchor_push(struct ANCHOR_V *v, struct ANCHOR_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct ANCHOR_T *)realloc(v->buffer, sizeof(struct ANCHOR_T) * v->cap);
	}
	v->buffer[v->size++] = t;
}

struct CHAIN_T
{
	b4i one;
	b4i more;
};

struct CHAIN_V
{
	u8i size;
	u8i cap;
	struct CHAIN_T *buffer;
};

void chain_push(struct CHAIN_V *v, struct CHAIN_T t)
{
	u8i ui;
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct CHAIN_T *)realloc(v->buffer, sizeof(struct CHAIN_T) * v->cap);
	}
	v->size++;
	for(ui = v->size - 1; ui >= 1; ui--)
	{
		v->buffer[ui] = v->buffer[ui - 1];
	}
	v->buffer[0] = t;
}

void initqrybin(BioSequence *seq, b4i ksize, b4i wsize, b4i bsize, b4i bi, struct MINIMIZER_V *minimizer_v, struct QRYBIN_V **qrybin_v)
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

float geta(struct ANCHOR_T anchor_j, struct ANCHOR_T anchor_i, b4i ksize)
{
	float tmp;
	tmp = ksize;
	if(tmp > anchor_i.y - anchor_j.y)
	{
		tmp = anchor_i.y - anchor_j.y;
	}
	if(tmp > anchor_i.x - anchor_j.x)
	{
		tmp = anchor_i.x - anchor_j.x;
	}
	return tmp;
}

float getb(struct ANCHOR_T anchor_j, struct ANCHOR_T anchor_i, b4i ksize)
{
	float tmp;
	tmp = (anchor_i.y - anchor_j.y) - (anchor_i.x - anchor_j.x);
	if(tmp != 0)
	{
		tmp = tmp < 0 ? 0 - tmp : tmp;
		tmp = 0.01 * ksize * tmp + 0.5 * (log(tmp) / log(2));
	}
	return tmp;
}

void updateqrybin(BioSequence *seq, b4i ksize, b4i wsize, b4i bsize, b4i bi, struct MINIMIZER_V *kmer_v, struct QRYBIN_V **bin_v)
{
	b4i pre_idx, min_idx, bj, bk, bl;
	u1i nucleobase[256];
	u8i key[4], mask[4], ui;
	struct MINIMIZER_T pre_elm[64], min_elm;
	struct QRYBIN_T qrybin_t;
	u8i uj;
	b4i maxp;
	float maxf, tmp;
	struct MINIMIZER_V *minimizer_v;
	struct ANCHOR_T anchor_t;
	struct ANCHOR_V *anchor_v;
	struct CHAIN_T chain_t;
	struct CHAIN_V *chain_v;
	b4i beg[seq->seq->size - ksize + 1], end[seq->seq->size - ksize + 1];
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
	minimizer_v = (struct MINIMIZER_V *)malloc(sizeof(struct MINIMIZER_V));
	minimizer_v->size = minimizer_v->cap = 0;
	minimizer_v->buffer = NULL;
	anchor_v = (struct ANCHOR_V *)malloc(sizeof(struct ANCHOR_V));
	anchor_v->size = anchor_v->cap = 0;
	anchor_v->buffer = NULL;
	chain_v = (struct CHAIN_V *)malloc(sizeof(struct CHAIN_V));
	chain_v->size = chain_v->cap = 0;
	chain_v->buffer = NULL;
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
				for(ui = 0; ui < kmer_v->size; ui++)
				{
					if(pre_elm[pre_idx].s == kmer_v->buffer[ui].s)
					{
						anchor_t.x = bj - ksize + 1;
						anchor_t.y = kmer_v->buffer[ui].ip & mask[3];
						anchor_t.f = ksize;
						anchor_t.p = -1;
						anchor_push(anchor_v, anchor_t);
					}
				}
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
							minimizer_push(minimizer_v, pre_elm[bl]);
						}
					}
					for(bl = 0; bl < pre_idx; bl++)
					{
						if(pre_elm[bl].s == min_elm.s)
						{
							minimizer_push(minimizer_v, pre_elm[bl]);
						}
					}
				}
				else
				{
					if(min_elm.s >= pre_elm[pre_idx].s)
					{
						minimizer_push(minimizer_v, pre_elm[pre_idx]);
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
									minimizer_push(minimizer_v, pre_elm[bl]);
								}
							}
							for(bl = 0; bl < pre_idx; bl++)
							{
								if(pre_elm[bl].s == min_elm.s)
								{
									minimizer_push(minimizer_v, pre_elm[bl]);
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
	maxf = ksize;
	maxp = 0;
	for(ui = 1; ui < anchor_v->size; ui++)
	{
		for(uj = 0; uj < ui; uj++)
		{
			if(anchor_v->buffer[uj].x >= anchor_v->buffer[ui].x || anchor_v->buffer[uj].y >= anchor_v->buffer[ui].y)
			{
				continue;
			}
			tmp = anchor_v->buffer[uj].f + geta(anchor_v->buffer[uj], anchor_v->buffer[ui], ksize) - getb(anchor_v->buffer[uj], anchor_v->buffer[ui], ksize);
			if(anchor_v->buffer[ui].f < tmp)
			{
				anchor_v->buffer[ui].f = tmp;
				anchor_v->buffer[ui].p = uj;
			}
		}
		if(maxf < anchor_v->buffer[ui].f)
		{
			maxf = anchor_v->buffer[ui].f;
			maxp = ui;
		}
	}
	chain_t.one = kmer_v->size - 1;
	chain_t.more = seq->seq->size - ksize;
	chain_push(chain_v, chain_t);
	while(maxp != -1)
	{
		chain_t.one = anchor_v->buffer[maxp].y;
		chain_t.more = anchor_v->buffer[maxp].x;
		chain_push(chain_v, chain_t);
		maxp = anchor_v->buffer[maxp].p;
	}
	chain_t.one = 0;
	chain_t.more = 0;
	chain_push(chain_v, chain_t);
	for(ui = 1; ui < chain_v->size; ui++)
	{
		if(ui != 1)
		{
			bj = chain_v->buffer[ui - 1].more;
			beg[bj] = chain_v->buffer[ui - 1].one;
			beg[bj] = beg[bj] - (beg[bj] % bsize);
			beg[bj] = beg[bj] / bsize;
			end[bj] = beg[bj];
		}
		for(bj = chain_v->buffer[ui - 1].more + 1; bj <= chain_v->buffer[ui].more - 1; bj++)
		{
			beg[bj] = chain_v->buffer[ui - 1].one;
			beg[bj] = beg[bj] - (beg[bj] % bsize);
			beg[bj] = beg[bj] / bsize;
			end[bj] = chain_v->buffer[ui].one;
			end[bj] = end[bj] - (end[bj] % bsize);
			end[bj] = end[bj] / bsize;
		}
		if(ui != chain_v->size - 1)
		{
			bj = chain_v->buffer[ui].more;
			end[bj] = chain_v->buffer[ui].one;
			end[bj] = end[bj] - (end[bj] % bsize);
			end[bj] = end[bj] / bsize;
			beg[bj] = end[bj];
		}
	}
	for(ui = 0; ui < minimizer_v->size; ui++)
	{
		qrybin_t.s = minimizer_v->buffer[ui].s;
		bj = minimizer_v->buffer[ui].ip & mask[3];
		for(bk = beg[bj]; bk <= end[bj]; bk++)
		{
			qrybin_push(bin_v[bk], qrybin_t);
		}
	}
	if(minimizer_v->size > 0)
	{
		free(minimizer_v->buffer);
	}
	free(minimizer_v);
	if(anchor_v->size > 0)
	{
		free(anchor_v->buffer);
	}
	free(anchor_v);
	if(chain_v->size > 0)
	{
		free(chain_v->buffer);
	}
	free(chain_v);
}

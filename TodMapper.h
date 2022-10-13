#include <math.h>
#include "bsalign/filereader.h"

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

struct MYKMER_T
{
	u8i s;
	u4i i;
	u4i p;
};

struct MYKMER_V
{
	u8i size;
	u8i cap;
	struct MYKMER_T *buffer;
	u4i *i;
	u4i *p;
	u4i bb[256];
	u4i bl[256];
};

struct MYKMER2_V
{
	u8i size;
	u8i cap;
	struct MYKMER_T *buffer;
};

void mykmer2_push(struct MYKMER2_V *v, struct MYKMER_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct MYKMER_T *)realloc(v->buffer, sizeof(struct MYKMER_T) * v->cap);
	}
	v->buffer[v->size++] = t;
}

void mykmer_push(struct MYKMER_V *v, struct MYKMER_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct MYKMER_T *)realloc(v->buffer, sizeof(struct MYKMER_T) * v->cap);
	}
	v->buffer[v->size++] = t;
}

void mykmer_push2(struct MYKMER_V *v, struct MYKMER_T t, u4i p)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct MYKMER_T *)realloc(v->buffer, sizeof(struct MYKMER_T) * v->cap);
	}
	v->buffer[v->size++] = t;
	t.p = p - t.p;
	mykmer_push(v, t);
}

struct MYANCHOR_T
{
	u8i x;
	u4i y;
	u4i bcnt;
};

struct MYANCHOR_V
{
	u8i size;
	u8i cap;
	struct MYANCHOR_T *buffer;
};

void myanchor_push(struct MYANCHOR_V *v, struct MYANCHOR_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct MYANCHOR_T *)realloc(v->buffer, sizeof(struct MYANCHOR_T) * v->cap);
	}
	v->buffer[v->size++] = t;
}

void myanchor_push2(struct MYANCHOR_V *v, struct MYANCHOR_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct MYANCHOR_T *)realloc(v->buffer, sizeof(struct MYANCHOR_T) * v->cap);
	}
	if(v->size == 0 || v->buffer[v->size - 1].bcnt <= t.bcnt)
	{
		v->buffer[v->size++] = t;
	}
}

void mystep1(const char *seqseqstring, string_size_t seqseqsize, b4i ksize, b4i wsize, b4i bi, struct MYKMER2_V *mykmer2_v, u1i *nucleobase, u8i *mask)
{
	b4i pre_idx, min_idx, bj, bk, bl;
	u1i z;
	u8i key[4];
	struct MYKMER_T pre_elm[64], min_elm;
	//
	pre_idx = min_idx = 0;
	key[0] = key[1] = 0;
	min_elm.s = UINT64_MAX;
	bk = 0;
	for(bj = 0; bj < 64; bj++)
	{
		pre_elm[bj].i = bi;
	}
	for(bj = 0; bj < seqseqsize; bj++)
	{
		nucleobase[128] = nucleobase[(int)seqseqstring[bj]];
		if(nucleobase[128] < 4)
		{
			key[0] = ((key[0] << 2) & mask[0]) | nucleobase[128];
			key[1] = ((key[1] >> 2) & mask[1]) | ((nucleobase[128] ^ 3ULL) << mask[2]);
			bk++;
			if(bk >= ksize)
			{
				z = key[0] <= key[1] ? 0 : 1;
				pre_elm[pre_idx].s = hash64(key[z], mask[0]);
				pre_elm[pre_idx].p = (bj - ksize + 1) >> mask[6];
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
							mykmer2_push(mykmer2_v, pre_elm[bl]);
						}
					}
					for(bl = 0; bl < pre_idx; bl++)
					{
						if(pre_elm[bl].s == min_elm.s)
						{
							mykmer2_push(mykmer2_v, pre_elm[bl]);
						}
					}
				}
				else
				{
					if(min_elm.s >= pre_elm[pre_idx].s)
					{
						mykmer2_push(mykmer2_v, pre_elm[pre_idx]);
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
									mykmer2_push(mykmer2_v, pre_elm[bl]);
								}
							}
							for(bl = 0; bl < pre_idx; bl++)
							{
								if(pre_elm[bl].s == min_elm.s)
								{
									mykmer2_push(mykmer2_v, pre_elm[bl]);
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

static void mystep12(const struct MYKMER_T *a, int n, struct MYKMER_V **mykmer_v)
{
	int i, mask = (1 << 14) - 1;
	for(i = 0; i < n; i++)
	{
		mykmer_push(mykmer_v[a[i].s & mask], a[i]);
	}
}

void mystep2(struct MYKMER_V *mykmer_vs, struct MYKMER_V **mykmer_v2, u8i *mask)
{
	u8i ui, uj, uk, ul;
	struct MYKMER_T mykmer_t;
	for(ui = mask[4] + 8; ui < 64; ui = ui + 8)
	{
		for(uj = 0; uj < mykmer_vs->size; uj++)
		{
			mykmer_push(mykmer_v2[(mykmer_vs->buffer[uj].s >> ui) & 255ULL], mykmer_vs->buffer[uj]);
		}
		mykmer_vs->size = 0;
		for(uj = 0; uj < 256; uj++)
		{
			for(uk = 0; uk < mykmer_v2[uj]->size; uk++)
			{
				mykmer_push(mykmer_vs, mykmer_v2[uj]->buffer[uk]);
			}
			mykmer_v2[uj]->size = 0;
		}
	}
	ui = mask[4];
	for(uj = 0; uj < mykmer_vs->size; uj++)
	{
		mykmer_push(mykmer_v2[(mykmer_vs->buffer[uj].s >> ui) & 255ULL], mykmer_vs->buffer[uj]);
	}
	mykmer_vs->size = 0;
	ul = 0;
	for(uj = 0; uj < 256; uj++)
	{
		mykmer_vs->bb[uj] = mykmer_vs->size;
		mykmer_vs->bl[uj] = 0;
		if(mykmer_v2[uj]->size)
		{
			uk = 0;
			mykmer_t.s = mykmer_v2[uj]->buffer[uk].s;
			mykmer_t.i = ul;
			mykmer_t.p = 0;
			mykmer_t.p++;
			mykmer_vs->i[ul] = mykmer_v2[uj]->buffer[uk].i;
			mykmer_vs->p[ul++] = mykmer_v2[uj]->buffer[uk].p;
			for(uk = 1; uk < mykmer_v2[uj]->size; uk++)
			{
				if(mykmer_t.s != mykmer_v2[uj]->buffer[uk].s)
				{
					mykmer_push(mykmer_vs, mykmer_t);
					mykmer_vs->bl[uj]++;
					mykmer_t.s = mykmer_v2[uj]->buffer[uk].s;
					mykmer_t.i = ul;
					mykmer_t.p = 0;
				}
				mykmer_t.p++;
				mykmer_vs->i[ul] = mykmer_v2[uj]->buffer[uk].i;
				mykmer_vs->p[ul++] = mykmer_v2[uj]->buffer[uk].p;
			}
			mykmer_push(mykmer_vs, mykmer_t);
			mykmer_vs->bl[uj]++;
			mykmer_v2[uj]->size = 0;
		}
	}
}

void mystep3(const char *seqseqstring, string_size_t seqseqsize, b4i ksize, b4i wsize, b4i bi, struct MYKMER_V **mykmer_v, u1i *nucleobase, u8i *mask)
{
        b4i pre_idx, min_idx, bj, bk, bl;
	u1i z;
        u8i key[4], p = (seqseqsize - ksize) >> mask[6];
        struct MYKMER_T pre_elm[64], min_elm;
	//
        pre_idx = min_idx = 0;
        key[0] = key[1] = 0;
        min_elm.s = UINT64_MAX;
        bk = 0;
	for(bj = 0; bj < 64; bj++)
	{
		pre_elm[bj].i = bi;
	}
        for(bj = 0; bj < seqseqsize; bj++)
        {
                nucleobase[128] = nucleobase[(int)seqseqstring[bj]];
                if(nucleobase[128] < 4)
                {
                        key[0] = ((key[0] << 2) & mask[0]) | nucleobase[128];
                        key[1] = ((key[1] >> 2) & mask[1]) | ((nucleobase[128] ^ 3ULL) << mask[2]);
                        bk++;
                        if(bk >= ksize)
                        {
				z = key[0] <= key[1] ? 0 : 1;
				pre_elm[pre_idx].s = hash64(key[z], mask[0]);
                                pre_elm[pre_idx].p = (bj - ksize + 1) >> mask[6];
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
                                                        mykmer_push2(mykmer_v[pre_elm[bl].s & 255ULL], pre_elm[bl], p);
                                                }
                                        }
                                        for(bl = 0; bl < pre_idx; bl++)
                                        {
                                                if(pre_elm[bl].s == min_elm.s)
                                                {
							mykmer_push2(mykmer_v[pre_elm[bl].s & 255ULL], pre_elm[bl], p);
                                                }
                                        }
                                }
                                else
                                {
                                        if(min_elm.s >= pre_elm[pre_idx].s)
                                        {
						mykmer_push2(mykmer_v[pre_elm[pre_idx].s & 255ULL], pre_elm[pre_idx], p);
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
									mykmer_push2(mykmer_v[pre_elm[bl].s & 255ULL], pre_elm[bl], p);
                                                                }
                                                        }
                                                        for(bl = 0; bl < pre_idx; bl++)
                                                        {
                                                                if(pre_elm[bl].s == min_elm.s)
                                                                {
									mykmer_push2(mykmer_v[pre_elm[bl].s & 255ULL], pre_elm[bl], p);
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

void mystep32(struct MYKMER_V *mykmer_v, b4i beg, b4i end)
{
        b4i left, right;
        struct MYKMER_T mykmer_t;
        left = beg;
        right = end;
        mykmer_t = mykmer_v->buffer[left];
        while(left != right)
        {
                for(right = right; right > left; right--)
                {
                        if((mykmer_v->buffer[right].s < mykmer_t.s) || (mykmer_v->buffer[right].s == mykmer_t.s && mykmer_v->buffer[right].p < mykmer_t.p))
                        {
                                mykmer_v->buffer[left] = mykmer_v->buffer[right];
                                break;
                        }
                }
                for(left = left; left < right; left++)
                {
                        if((mykmer_v->buffer[left].s > mykmer_t.s) || (mykmer_v->buffer[left].s == mykmer_t.s && mykmer_v->buffer[left].p > mykmer_t.p))
                        {
                                mykmer_v->buffer[right] = mykmer_v->buffer[left];
                                break;
                        }
                }
        }
        mykmer_v->buffer[left] = mykmer_t;
        if(left - beg > 1)
        {
                mystep32(mykmer_v, beg, left - 1);
        }
        if(end - right > 1)
        {
                mystep32(mykmer_v, right + 1, end);
        }
}

void mystep33(struct MYKMER_V **mykmer_v)
{
	u8i ui, uj, uk, ul;
	struct MYKMER_T mykmer_t;
	for(ui = 0; ui < 256; ui++)
	{
		for(uj = 0; uj < mykmer_v[ui]->size; uj++)
		{
			uk = 256 + ((mykmer_v[ui]->buffer[uj].s >> 8) & 255ULL);
			mykmer_push(mykmer_v[uk], mykmer_v[ui]->buffer[uj]);
		}
		mykmer_v[ui]->size = 0;
		for(uk = 256; uk < 512; uk++)
		{
			if(mykmer_v[uk]->size > 1)
			{
				mystep32(mykmer_v[uk], 0, mykmer_v[uk]->size - 1);
				mykmer_t.s = mykmer_v[uk]->buffer[0].s;
				mykmer_t.i = mykmer_v[uk]->buffer[0].i;
				mykmer_t.p = mykmer_v[uk]->buffer[0].p;
				for(ul = 1; ul < mykmer_v[uk]->size; ul++)
				{
					if(mykmer_t.s != mykmer_v[uk]->buffer[ul].s || mykmer_t.p != mykmer_v[uk]->buffer[ul].p)
					{
						if(!mykmer_t.i)
						{
							mykmer_push(mykmer_v[0], mykmer_t);
						}
						mykmer_t.s = mykmer_v[uk]->buffer[ul].s;
						mykmer_t.i = mykmer_v[uk]->buffer[ul].i;
						mykmer_t.p = mykmer_v[uk]->buffer[ul].p;
					}
					else
					{
						mykmer_t.i = mykmer_t.i == mykmer_v[uk]->buffer[ul].i ? mykmer_t.i : 0;
					}
				}
				if(!mykmer_t.i)
				{
					mykmer_push(mykmer_v[0], mykmer_t);
				}
			}
			mykmer_v[uk]->size = 0;
		}
	}
}

b4i mystep34(u8i s, struct MYKMER_V *mykmer_v, u4i bb, u4i bl)
{
	b4i low, mid, high, bi;
	low = bb;
	mid = low + bl;
	high = mid - 1;
	while(low <= high)
	{
		bi = low + high;
		mid = bi % 2 == 0 ? bi / 2 : (bi - 1) / 2;
		if(mykmer_v->buffer[mid].s < s)
		{
			low = mid + 1;
		}
		else if(mykmer_v->buffer[mid].s > s)
		{
			high = mid - 1;
		}
		else
		{
			break;
		}
	}
	if(low > high)
	{
		mid = bb + bl;
	}
	return mid;
}

void mystep4(struct MYKMER_V *mykmer_v2, b4i kmax, struct MYANCHOR_V **myanchor_v, struct MYKMER_V **mykmer_vs, u8i *mask)
{
	u8i ui, uj, uk, ul, um, un;
	struct MYANCHOR_T myanchor_t;
	for(ui = 0; ui < mykmer_v2->size; ui++)
	{
		uj = mykmer_v2->buffer[ui].s & (mask[5] - 1);
		uk = (mykmer_v2->buffer[ui].s >> mask[4]) & 255ULL;
		ul = mykmer_vs[uj]->bl[uk] ? (u8i)mystep34(mykmer_v2->buffer[ui].s, mykmer_vs[uj], mykmer_vs[uj]->bb[uk], mykmer_vs[uj]->bl[uk]) : mykmer_vs[uj]->bb[uk] + mykmer_vs[uj]->bl[uk];
		if(ul != mykmer_vs[uj]->bb[uk] + mykmer_vs[uj]->bl[uk] && mykmer_vs[uj]->buffer[ul].p <= (u4i)kmax)
		{
			myanchor_t.y = mykmer_v2->buffer[ui].p;
			myanchor_t.bcnt = 0;
			um = mykmer_vs[uj]->buffer[ul].i + mykmer_vs[uj]->buffer[ul].p;
			for(un = mykmer_vs[uj]->buffer[ul].i; un < um; un++)
			{
				myanchor_t.x = mykmer_vs[uj]->p[un] <= myanchor_t.y ? 0 : ((u8i)mykmer_vs[uj]->i[un] << 32) | (mykmer_vs[uj]->p[un] - myanchor_t.y);
				myanchor_push(myanchor_v[myanchor_t.x & 255ULL], myanchor_t);
			}
		}
	}
	mykmer_v2->size = 0;
}

void mystep45(struct MYANCHOR_V *myanchor_v, b4i beg, b4i end)
{
        b4i left, right;
        struct MYANCHOR_T myanchor_t;
        left = beg;
        right = end;
        myanchor_t = myanchor_v->buffer[left];
        while(left != right)
        {
                for(right = right; right > left; right--)
                {
			if(myanchor_v->buffer[right].x < myanchor_t.x || (myanchor_v->buffer[right].x == myanchor_t.x && myanchor_v->buffer[right].y < myanchor_t.y))
                        {
                                myanchor_v->buffer[left] = myanchor_v->buffer[right];
                                break;
                        }
                }
                for(left = left; left < right; left++)
                {
			if(myanchor_v->buffer[left].x > myanchor_t.x || (myanchor_v->buffer[left].x == myanchor_t.x && myanchor_v->buffer[left].y > myanchor_t.y))
                        {
                                myanchor_v->buffer[right] = myanchor_v->buffer[left];
                                break;
                        }
                }
        }
        myanchor_v->buffer[left] = myanchor_t;
        if(left - beg > 1)
        {
                mystep45(myanchor_v, beg, left - 1);
        }
        if(end - right > 1)
        {
                mystep45(myanchor_v, right + 1, end);
        }
}

void mystep5(struct MYANCHOR_V **myanchor_vs, b4i bsize, b4i bmin, b4i mmin, FILE *outfp, struct MYANCHOR_V **myanchor_v2, BioSequence *seq, u8i *mask)
{
	b4i bi;
	u4i mcnt;
	u8i ui, uj, uk;
	struct MYANCHOR_T myanchor_t;
	for(ui = 0; ui < 256; ui++)
	{
		for(uj = 0; uj < myanchor_vs[ui]->size; uj++)
		{
			myanchor_push(myanchor_v2[(myanchor_vs[ui]->buffer[uj].x >> 8) & 255ULL], myanchor_vs[ui]->buffer[uj]);
		}
		myanchor_vs[ui]->size = 0;
		for(uj = 0; uj < 256; uj++)
		{
			if(myanchor_v2[uj]->size > (u8i)bmin * (u8i)mmin)
			{
				mystep45(myanchor_v2[uj], 0, myanchor_v2[uj]->size - 1);
				for(uk = 0; uk < myanchor_v2[uj]->size; uk++)
				{
					myanchor_push(myanchor_vs[0], myanchor_v2[uj]->buffer[uk]);
				}
			}
			myanchor_v2[uj]->size = 0;
		}
	}
	ui = 0;
	myanchor_t.x = myanchor_vs[0]->buffer[ui].x;
	myanchor_t.y = myanchor_vs[0]->buffer[ui].y;
	myanchor_t.bcnt = 0;
	mcnt = 1;
	for(ui = 1; ui < myanchor_vs[0]->size; ui++)
	{
		if(myanchor_t.x == myanchor_vs[0]->buffer[ui].x)
		{
			if(myanchor_t.y == myanchor_vs[0]->buffer[ui].y)
			{
				mcnt++;
			}
			else
			{
				if(mcnt >= (u4i)mmin)
				{
					myanchor_t.bcnt++;
				}
				myanchor_t.y = myanchor_vs[0]->buffer[ui].y;
				mcnt = 1;
			}
		}
		else
		{
			if(mcnt >= (u4i)mmin)
			{
				myanchor_t.bcnt++;
			}
			if(myanchor_t.bcnt >= (u4i)bmin)
			{
				myanchor_push2(myanchor_v2[0], myanchor_t);
			}
			myanchor_t.x = myanchor_vs[0]->buffer[ui].x;
			myanchor_t.y = myanchor_vs[0]->buffer[ui].y;
			myanchor_t.bcnt = 0;
			mcnt = 1;
		}
	}
	if(mcnt >= (u4i)mmin)
	{
		myanchor_t.bcnt++;
	}
	if(myanchor_t.bcnt >= (u4i)bmin)
	{
		myanchor_push2(myanchor_v2[0], myanchor_t);
	}
	myanchor_vs[0]->size = 0;
	if(myanchor_v2[0]->size > 0)
	{
		bi = myanchor_v2[0]->size - 1;
		fprintf(outfp, "%s %llu %llu %u %u\n", seq->tag->string, (myanchor_v2[0]->buffer[bi].x >> 32) & mask[3], (myanchor_v2[0]->buffer[bi].x & mask[3]) * bsize, (myanchor_v2[0]->buffer[bi].y + 1) * bsize, myanchor_v2[0]->buffer[bi].bcnt);
		for(bi = bi - 1; bi >= 0 && myanchor_v2[0]->buffer[bi].bcnt == myanchor_v2[0]->buffer[bi + 1].bcnt; bi--)
		{
			fprintf(outfp, "%s %llu %llu %u %u\n", seq->tag->string, (myanchor_v2[0]->buffer[bi].x >> 32) & mask[3], (myanchor_v2[0]->buffer[bi].x & mask[3]) * bsize, (myanchor_v2[0]->buffer[bi].y + 1) * bsize, myanchor_v2[0]->buffer[bi].bcnt);
		}
		myanchor_v2[0]->size = 0;
	}
}

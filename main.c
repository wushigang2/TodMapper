#include "TodMapper.h"

int usage_cnt()
{
	return 1;
}

int main_cnt()
{
	if(_DEBUG_LOG_)
	{
	}
	return 0;
}

int usage_map()
{
	fprintf(stdout,
	"usage: TodMapper map [options] -q <input file> [-Q <input file>] -r <input file> -o <output file>\n"
	" -q <string> multi sequence file, [NULL]\n"
	" -r <string> genome file, [NULL]\n"
	" -Q <string> the second multi sequence file, [NULL]\n"
	" -o <string> output file, [NULL]\n"
	" -k <int>    kmer size, [15]\n"
	" -w <int>    window size, [10]\n"
	" -b <int>    bin size, [256]\n"
	" -K <int>    filter high frequency minimizers, maybe repetitive, [4096]\n"
	" -e <int>    min minimizer depth of a valid edge, [2]\n"
	" -E <int>    top minimizer depth of a valid edge, [128]\n"
	" -N <int>    retain at most N secondary alignments, [9]\n"
	" -d <int>    depth of the multisequence, [10]\n"
	" -g <int>    gaps of a bin alignment, [5]\n"
	" -v          verbose\n"
	"\n"
	"for example: map multi sequence qry.fa to genome ref.fa\n"
	"TodMapper map -q qry.fa -r ref.fa -o output.txt\n"
	"for example: calculate overlap of multi sequence qry.fa and the second multi sequence qry2.fa on genome ref.fa\n"
	"TodMapper map -q qry.fa -Q qry2.fa -r ref.fa -o output.txt\n"
	);
	return 1;
}

int main_map(int argc, char **argv)
{
	FileReader *fr, *fr2;
	BioSequence *seq, *seq2;
	char **qryfn, **reffn, **qryfn2;
	FILE *outfp;
	char *outfn;
	b4i c, ksize, wsize, bsize, kmax, bmin, mmin, mine, tope, outn, depth, bgap, ver, bi, bj;
	u1i nucleobase[256];
	u8i mask[8], ui;
	struct MYKMER_V **mykmer_vs;
	struct MYKMER_V **mykmer_v2;
	struct MYKMER2_V *mykmer2_v;
	struct MYANCHOR_V **myanchor_vs;
	struct MYANCHOR_V **myanchor_v2;
	struct TMP0_V *tmp0_v;
	struct TMP1_V *tmp1_v;
	struct TMP2_V *tmp2_v, *tmp2_a, *tmp2_b;
	struct TMP3_V *tmp3_v;
	struct TMP4_V *tmp4_v, *tmp4_a, *tmp4_b;
	struct TMP8_V *tmp8_v;
	//
	qryfn = (char **)malloc(sizeof(char *));
	qryfn[0] = NULL;
	reffn = (char **)malloc(sizeof(char *));
	reffn[0] = NULL;
	qryfn2 = (char **)malloc(sizeof(char *));
	qryfn2[0] = NULL;
	outfn = NULL;
	ksize = 15;
	wsize = 10;
	bsize = 256;
	kmax = 4096;
	bmin = 6;
	mmin = 1;
	mine = 2;
	tope = 128;
	outn = 9;
	depth = 10;
	bgap = 5;
	ver = 0;
	while((c = getopt(argc, argv, "hq:r:Q:o:k:w:b:K:B:M:e:E:N:d:g:v")) != -1)
	{
		switch(c)
		{
			case 'q':
				qryfn[0] = optarg;
				break;
			case 'r':
				reffn[0] = optarg;
				break;
			case 'Q':
				qryfn2[0] = optarg;
				break;
			case 'o':
				outfn = optarg;
				break;
			case 'k':
				ksize = atoi(optarg);
				break;
			case 'w':
				wsize = atoi(optarg);
				break;
			case 'b':
				bsize = atoi(optarg);
				break;
			case 'K':
				kmax = atoi(optarg);
				break;
			case 'B':
				bmin = atoi(optarg);
				break;
			case 'M':
				mmin = atoi(optarg);
				break;
			case 'e':
				mine = atoi(optarg);
				break;
			case 'E':
				tope = atoi(optarg);
				break;
			case 'N':
				outn = atoi(optarg);
				break;
			case 'd':
				depth = atoi(optarg);
				break;
			case 'g':
				bgap = atoi(optarg);
				break;
			case 'v':
				ver++;
				break;
			default:
				return usage_map();
		}
	}
	if((optind == argc && optind == 1) || (optind != argc))
	{
		return usage_map();
	}
	//
	for(bi = 0; bi < 128; bi++)
	{
		nucleobase[bi] = 4;
	}
	nucleobase[65] = nucleobase[97] = 0;
	nucleobase[67] = nucleobase[99] = 1;
	nucleobase[71] = nucleobase[103] = 2;
	nucleobase[84] = nucleobase[85] = nucleobase[116] = nucleobase[117] = 3;
	mask[0] = (1ULL << (ksize * 2)) - 1;
	mask[1] = (1ULL << 62) - 1;
	mask[2] = (ksize - 1) * 2;
	mask[3] = (1ULL << 32) - 1;
	mask[4] = 14;
	mask[5] = 1ULL << mask[4];
	mask[6] = (u8i)(log(bsize) / log(2));
	mykmer_vs = (struct MYKMER_V **)calloc(mask[5], sizeof(struct MYKMER_V *));
	for(ui = 0; ui < mask[5]; ui++)
	{
		mykmer_vs[ui] = (struct MYKMER_V *)malloc(sizeof(struct MYKMER_V));
		mykmer_vs[ui]->size = mykmer_vs[ui]->cap = 0;
		mykmer_vs[ui]->buffer = NULL;
		mykmer_vs[ui]->i = NULL;
		mykmer_vs[ui]->p = NULL;
	}
	mykmer_v2 = (struct MYKMER_V **)calloc(512, sizeof(struct MYKMER_V *));
	for(ui = 0; ui < 512; ui++)
	{
		mykmer_v2[ui] = (struct MYKMER_V *)malloc(sizeof(struct MYKMER_V));
		mykmer_v2[ui]->size = mykmer_v2[ui]->cap = 0;
		mykmer_v2[ui]->buffer = NULL;
		mykmer_v2[ui]->i = NULL;
		mykmer_v2[ui]->p = NULL;
	}
	mykmer2_v = (struct MYKMER2_V *)malloc(sizeof(struct MYKMER2_V));
	mykmer2_v->size = mykmer2_v->cap = 0;
	mykmer2_v->buffer = NULL;
	myanchor_vs = (struct MYANCHOR_V **)calloc(256, sizeof(struct MYANCHOR_V *));
	myanchor_v2 = (struct MYANCHOR_V **)calloc(256, sizeof(struct MYANCHOR_V *));
	for(ui = 0; ui < 256; ui++)
	{
		myanchor_vs[ui] = (struct MYANCHOR_V *)malloc(sizeof(struct MYANCHOR_V));
		myanchor_vs[ui]->size = myanchor_vs[ui]->cap = 0;
		myanchor_vs[ui]->buffer = NULL;
		myanchor_v2[ui] = (struct MYANCHOR_V *)malloc(sizeof(struct MYANCHOR_V));
		myanchor_v2[ui]->size = myanchor_v2[ui]->cap = 0;
		myanchor_v2[ui]->buffer = NULL;
	}
	tmp0_v = (struct TMP0_V *)malloc(sizeof(struct TMP0_V));
	tmp1_v = (struct TMP1_V *)malloc(sizeof(struct TMP1_V));
	tmp2_a = (struct TMP2_V *)malloc(sizeof(struct TMP2_V));
	tmp2_b = (struct TMP2_V *)malloc(sizeof(struct TMP2_V));
	tmp2_v = (struct TMP2_V *)malloc(sizeof(struct TMP2_V));
	tmp3_v = (struct TMP3_V *)malloc(sizeof(struct TMP3_V));
	tmp4_v = (struct TMP4_V *)malloc(sizeof(struct TMP4_V));
	tmp4_a = (struct TMP4_V *)malloc(sizeof(struct TMP4_V));
	tmp4_b = (struct TMP4_V *)malloc(sizeof(struct TMP4_V));
	tmp8_v = (struct TMP8_V *)malloc(sizeof(struct TMP8_V));
	tmp0_v->size = tmp1_v->size = tmp2_a->size = tmp2_b->size = tmp2_v->size = tmp3_v->size = tmp4_v->size = tmp4_a->size = tmp4_b->size = tmp8_v->size = 0;
	tmp0_v->cap = tmp1_v->cap = tmp2_a->cap = tmp2_b->cap = tmp2_v->cap = tmp3_v->cap = tmp4_v->cap = tmp4_a->cap = tmp4_b->cap = tmp8_v->cap = 0;
	//
	clock_t start0, end0, sme0 = 0, start01, end01, sme01 = 0, start02, end02, sme02 = 0;
	//clock_t start0, end0, sme0 = 0;
	double usersys0, usersys01, usersys02;
	//double usersys0;
	start0 = clock();
	fr = open_all_filereader(1, reffn, 0);
	free(reffn);
	seq = init_biosequence();
	start01 = clock();
	bi = 0;
	while(readseq_filereader(fr, seq))
	{
		mystep1(seq->seq->string, seq->seq->size, ksize, wsize, bi, mykmer2_v, nucleobase, mask);
		mystep12(mykmer2_v->buffer, mykmer2_v->size, mykmer_vs);
		mykmer2_v->size = 0;
		tmp0_push(tmp0_v, seq->tag->string);
		bi++;
	}
	if(mykmer2_v->buffer)
	{
		free(mykmer2_v->buffer);
	}
	free(mykmer2_v);
	end01 = clock();
	sme01 += end01 - start01;
	free_biosequence(seq);
	close_filereader(fr);
	start02 = clock();
	for(ui = 0; ui < mask[5]; ui++)
	{
		mykmer_vs[ui]->i = (u4i *)calloc(mykmer_vs[ui]->size, sizeof(u4i));
		mykmer_vs[ui]->p = (u4i *)calloc(mykmer_vs[ui]->size, sizeof(u4i));
		mystep2(mykmer_vs[ui], mykmer_v2, mask);
	}
	end02 = clock();
	sme02 += end02 - start02;
	end0 = clock();
	sme0 += end0 - start0;
	usersys0 = (double)sme0 / CLOCKS_PER_SEC;
	usersys01 = (double)sme01 / CLOCKS_PER_SEC;
	usersys02 = (double)sme02 / CLOCKS_PER_SEC;
	fprintf(stderr, "step1 = %f sec.\n", usersys01);
	fprintf(stderr, "step2 = %f sec.\n", usersys02);
	fprintf(stderr, "total = %f sec.\n", usersys0);
	//
	outfp = fopen(outfn, "w");
	if(qryfn2[0])
	{
		fr = open_all_filereader(1, qryfn, 0);
		fr2 = open_all_filereader(1, qryfn2, 0);
		free(qryfn);
		free(qryfn2);
		seq = init_biosequence();
		seq2 = init_biosequence();
		while(readseq_filereader(fr, seq))
		{
			get_tmp1(seq->seq->string, seq->seq->size, ksize, wsize, tmp1_v, nucleobase, mask);
			tmp1_to_tmp2_one(tmp1_v, tmp2_v);
			for(bi = 1; bi < depth; bi++)
			{
				readseq_filereader(fr, seq);
				get_tmp1(seq->seq->string, seq->seq->size, ksize, wsize, tmp1_v, nucleobase, mask);
				tmp1_to_tmp2_one(tmp1_v, tmp2_v);
			}
			if(tmp2_v->size > 1)
			{
				sort_tmp2_x_y(tmp2_v, 0, tmp2_v->size - 1);
			}
			tmp2_to_tmp3(tmp2_v, mine - 1, tmp3_v);
			if(tmp3_v->size > 1)
			{
				sort_tmp3_cnt(tmp3_v, 0, tmp3_v->size - 1);
			}
			tmp3_to_tmp2(tmp3_v, kmax, tope, tmp2_v, tmp2_a, tmp2_b, mykmer_vs, mask);
			if(tmp2_v->size > 1)
			{
				sort_tmp2_x(tmp2_v, 0, tmp2_v->size - 1);
			}
			tmp4_a->size = 0;
			tmp2_to_tmp4(tmp2_v, bgap, tmp4_a, tmp1_v, mask[3]);
			if(tmp4_a->size > 1)
			{
				sort_tmp4_kcnt_blen(tmp4_a, 0, tmp4_a->size - 1);
			}
			//
			readseq_filereader(fr2, seq2);
			get_tmp1(seq2->seq->string, seq2->seq->size, ksize, wsize, tmp1_v, nucleobase, mask);
			tmp1_to_tmp2_one(tmp1_v, tmp2_v);
			for(bi = 1; bi < depth; bi++)
			{
				readseq_filereader(fr2, seq2);
				get_tmp1(seq2->seq->string, seq2->seq->size, ksize, wsize, tmp1_v, nucleobase, mask);
				tmp1_to_tmp2_one(tmp1_v, tmp2_v);
			}
			if(tmp2_v->size > 1)
			{
				sort_tmp2_x_y(tmp2_v, 0, tmp2_v->size - 1);
			}
			tmp2_to_tmp3(tmp2_v, mine - 1, tmp3_v);
			if(tmp3_v->size > 1)
			{
				sort_tmp3_cnt(tmp3_v, 0, tmp3_v->size - 1);
			}
			tmp3_to_tmp2(tmp3_v, kmax, tope, tmp2_v, tmp2_a, tmp2_b, mykmer_vs, mask);
			if(tmp2_v->size > 1)
			{
				sort_tmp2_x(tmp2_v, 0, tmp2_v->size - 1);
			}
			tmp4_b->size = 0;
			tmp2_to_tmp4(tmp2_v, bgap, tmp4_b, tmp1_v, mask[3]);
			if(tmp4_b->size > 1)
			{
				sort_tmp4_kcnt_blen(tmp4_b, 0, tmp4_b->size - 1);
			}
			//
			tmp4_to_tmp8(tmp4_a, tmp4_b, outn + 1, tmp8_v, mask[6]);
			if(tmp8_v->size > 1)
			{
				sort_tmp8_kcnt_blen(tmp8_v, 0, tmp8_v->size - 1);
			}
			put_tmp8(tmp8_v, outfp, seq->tag->string, seq2->tag->string, tmp0_v);
		}
		free_biosequence(seq);
		free_biosequence(seq2);
		close_filereader(fr);
		close_filereader(fr2);
	}
	else
	{
	clock_t start, end, sme = 0, start1, end1, sme1 = 0, start2, end2, sme2 = 0, start3, end3, sme3 = 0, start4, end4, sme4 = 0;
	//clock_t start, end, sme = 0;
	double usersys, usersys1, usersys2, usersys3, usersys4;
	//double usersys;
	start = clock();
	fr = open_all_filereader(1, qryfn, 0);
	free(qryfn);
	seq = init_biosequence();
	start1 = clock();
	while(readseq_filereader(fr, seq))
	{
		//mystep3(seq->seq->string, seq->seq->size, ksize, wsize, 1, mykmer_v2, nucleobase, mask);
		bj = seq->seq->size;
		get_tmp1(seq->seq->string, seq->seq->size, ksize, wsize, tmp1_v, nucleobase, mask);
		tmp1_to_tmp2_one(tmp1_v, tmp2_v);
		for(bi = 1; bi < depth; bi++)
		{
			readseq_filereader(fr, seq);
			//mystep3(seq->seq->string, seq->seq->size, ksize, wsize, bi + 1, mykmer_v2, nucleobase, mask);
			bj = bj >= seq->seq->size ? bj : seq->seq->size;
			get_tmp1(seq->seq->string, seq->seq->size, ksize, wsize, tmp1_v, nucleobase, mask);
			tmp1_to_tmp2_one(tmp1_v, tmp2_v);
		}
		if(tmp2_v->size > 1)
		{
			sort_tmp2_x_y(tmp2_v, 0, tmp2_v->size - 1);
		}
		/*fprintf(stdout, "sort_tmp2_x_y\n");
		for(ui = 0; ui < tmp2_v->size; ui++)
		{
			fprintf(stdout, "%llu x=%llu y=%llu\n", ui, tmp2_v->buffer[ui].x, tmp2_v->buffer[ui].y);
		}*/
		end1 = clock();
		sme1 += end1 - start1;
		start2 = clock();
		//mystep33(mykmer_v2);
		tmp2_to_tmp3(tmp2_v, mine - 1, tmp3_v);
		if(tmp3_v->size > 1)
		{
			sort_tmp3_cnt(tmp3_v, 0, tmp3_v->size - 1);
		}
		/*fprintf(stdout, "sort_tmp3_cnt\n");
		for(ui = 0; ui < tmp3_v->size; ui++)
		{
			fprintf(stdout, "%llu x=%llu y=%llu cnt=%llu\n", ui, tmp3_v->buffer[ui].x, tmp3_v->buffer[ui].y, tmp3_v->buffer[ui].cnt);
		}*/
		end2 = clock();
		sme2 += end2 - start2;
		start3 = clock();
		//mystep4(mykmer_v2[0], kmax, myanchor_vs, mykmer_vs, mask);
		tmp3_to_tmp2(tmp3_v, kmax, tope, tmp2_v, tmp2_a, tmp2_b, mykmer_vs, mask);
		if(tmp2_v->size > 1)
		{
			sort_tmp2_x(tmp2_v, 0, tmp2_v->size - 1);
		}
		/*fprintf(stdout, "sort_tmp1_s\n");
		for(ui = 0; ui < tmp1_v->size; ui++)
		{
			fprintf(stdout, "%llu ip=%llu\n", ui, tmp1_v->buffer[ui].s);
		}*/
		end3 = clock();
		sme3 += end3 - start3;
		start4 = clock();
		//mystep5(myanchor_vs, bsize, bmin, mmin, outfp, myanchor_v2, seq, mask);
		/*tmp1_to_tmp2_two(tmp1_v, tmp2_v, bj >> mask[6]);
		if(tmp2_v->size > 1)
		{
			sort_tmp2_y_x(tmp2_v, 0, tmp2_v->size - 1);
		}*/
		tmp2_to_tmp4(tmp2_v, bgap, tmp4_v, tmp1_v, mask[3]);
		if(tmp4_v->size > 1)
		{
			sort_tmp4_kcnt_blen(tmp4_v, 0, tmp4_v->size - 1);
		}
		/*fprintf(stdout, "sort_tmp2_x\n");
		for(ui = 0; ui < tmp2_v->size; ui++)
		{
			fprintf(stdout, "%llu ip=%llu score=%llu\n", ui, tmp2_v->buffer[ui].x, tmp2_v->buffer[ui].y);
		}*/
		/*if(tmp2_v->size > 0)
		{
			put_tmp2(tmp2_v, outn + 1, outfp, seq->tag->string, tmp0_v, (bj >> mask[6]) + 1, mask);
		}*/
		if(tmp4_v->size > 0)
		{
			put_tmp4(tmp4_v, outn + 1, outfp, seq->tag->string, tmp0_v, mask[6]);
		}
		end4 = clock();
		sme4 += end4 - start4;
		start1 = clock();
	}
	end1 = clock();
	sme1 += end1 - start1;
	free_biosequence(seq);
	close_filereader(fr);
	end = clock();
	sme += end - start;
	usersys1 = (double)sme1 / CLOCKS_PER_SEC;
	usersys2 = (double)sme2 / CLOCKS_PER_SEC;
	usersys3 = (double)sme3 / CLOCKS_PER_SEC;
	usersys4 = (double)sme4 / CLOCKS_PER_SEC;
	usersys = (double)sme / CLOCKS_PER_SEC;
	fprintf(stderr, "step1 = %f sec.\n", usersys1);
	fprintf(stderr, "step2 = %f sec.\n", usersys2);
	fprintf(stderr, "step3 = %f sec.\n", usersys3);
	fprintf(stderr, "step4 = %f sec.\n", usersys4);
	fprintf(stderr, "total = %f sec.\n", usersys);
	}
	fclose(outfp);
	//
	for(ui = 0; ui < mask[5]; ui++)
	{
		if(mykmer_vs[ui]->buffer)
		{
			free(mykmer_vs[ui]->buffer);
		}
		if(mykmer_vs[ui]->i)
		{
			free(mykmer_vs[ui]->i);
		}
		if(mykmer_vs[ui]->p)
		{
			free(mykmer_vs[ui]->p);
		}
		free(mykmer_vs[ui]);
	}
	free(mykmer_vs);
	for(ui = 0; ui < 512; ui++)
	{
		if(mykmer_v2[ui]->buffer)
		{
			free(mykmer_v2[ui]->buffer);
		}
		if(mykmer_v2[ui]->i)
		{
			free(mykmer_v2[ui]->i);
		}
		if(mykmer_v2[ui]->p)
		{
			free(mykmer_v2[ui]->p);
		}
		free(mykmer_v2[ui]);
	}
	free(mykmer_v2);
	for(ui = 0; ui < 256; ui++)
	{
		if(myanchor_vs[ui]->buffer)
		{
			free(myanchor_vs[ui]->buffer);
		}
		free(myanchor_vs[ui]);
		if(myanchor_v2[ui]->buffer)
		{
			free(myanchor_v2[ui]->buffer);
		}
		free(myanchor_v2[ui]);
	}
	free(myanchor_vs);
	free(myanchor_v2);
	if(tmp0_v->buffer)
	{
		for(ui = 0; ui < tmp0_v->size; ui++)
		{
			free(tmp0_v->buffer[ui].name);
		}
		free(tmp0_v->buffer);
	}
	free(tmp0_v);
	if(tmp1_v->buffer)
	{
		free(tmp1_v->buffer);
	}
	free(tmp1_v);
	if(tmp2_a->buffer)
	{
		free(tmp2_a->buffer);
	}
	free(tmp2_a);
	if(tmp2_b->buffer)
	{
		free(tmp2_b->buffer);
	}
	free(tmp2_b);
	if(tmp2_v->buffer)
	{
		free(tmp2_v->buffer);
	}
	free(tmp2_v);
	if(tmp3_v->buffer)
	{
		free(tmp3_v->buffer);
	}
	free(tmp3_v);
	if(tmp4_v->buffer)
	{
		free(tmp4_v->buffer);
	}
	free(tmp4_v);
	if(tmp4_a->buffer)
	{
		free(tmp4_a->buffer);
	}
	free(tmp4_a);
	if(tmp4_b->buffer)
	{
		free(tmp4_b->buffer);
	}
	free(tmp4_b);
	if(tmp8_v->buffer)
	{
		free(tmp8_v->buffer);
	}
	free(tmp8_v);
	//
	return 0;
}

int usage()
{
	fprintf(stdout,
	"program: TodMapper\n"
	"version: %s\n"
	"author : Shigang Wu <wushigang@caas.cn>\n"
	"usage  : TodMapper <cmd> [options]\n"
	"\n"
	"commands:\n"
	" map      map a file\n"
	, TOSTR(VERSION)
	);
	return 1;
}

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		return usage();
	}
	if(strcasecmp("cnt", argv[1]) == 0) return main_cnt(argc - 1, argv + 1);
	if(strcasecmp("map", argv[1]) == 0) return main_map(argc - 1, argv + 1);
	if(strcasecmp("-h", argv[1]) == 0) return usage();
	if(strcasecmp("--help", argv[1]) == 0) return usage();
	fprintf(stderr, " -- unknown command '%s' -- \n", argv[1]);
	return 1;
}

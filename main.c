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
	"usage: TodMapper map [options] <input file>\n"
	" -q <string> multi sequence file, [NULL]\n"
	" -r <string> genome file, [NULL]\n"
	" -o <string> output file, [NULL]\n"
	" -k <int>    kmer size, [15]\n"
	" -w <int>    window size, [10]\n"
	" -b <int>    bin size, [256]\n"
	" -K <int>    filter high frequency minimizers, maybe repetitive, [256]\n"
	" -B <int>    min bins of a match read, [4]\n"
	" -M <int>    min minimizers of a match bin, [2]\n"
	" -v          verbose\n"
	"\n"
	"for example: map multi sequence qry.fa to genome ref.fa\n"
	"TodMapper map -q qry.fa -r ref.fa -o output.txt\n"
	);
	return 1;
}

int main_map(int argc, char **argv)
{
	FileReader *fr;
	BioSequence *seq;
	char **qryfn, **reffn;
	FILE *outfp;
	char *outfn;
	b4i c, ksize, wsize, bsize, kmax, bmin, mmin, ver, bi;
	u1i nucleobase[256];
	u8i mask[8], ui;
	struct MYKMER_V **mykmer_vs;
	struct MYKMER_V **mykmer_v2;
	struct MYANCHOR_V **myanchor_vs;
	struct MYANCHOR_V **myanchor_v2;
	//
	qryfn = (char **)malloc(sizeof(char *));
	qryfn[0] = NULL;
	reffn = (char **)malloc(sizeof(char *));
	reffn[0] = NULL;
	outfn = NULL;
	ksize = 15;
	wsize = 10;
	bsize = 256;
	kmax = 256;
	bmin = 4;
	mmin = 2;
	ver = 0;
	while((c = getopt(argc, argv, "hq:r:o:k:w:b:K:B:M:v")) != -1)
	{
		switch(c)
		{
			case 'q':
				qryfn[0] = optarg;
				break;
			case 'r':
				reffn[0] = optarg;
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
	mykmer_v2 = (struct MYKMER_V **)calloc(256, sizeof(struct MYKMER_V *));
	myanchor_vs = (struct MYANCHOR_V **)calloc(256, sizeof(struct MYANCHOR_V *));
	myanchor_v2 = (struct MYANCHOR_V **)calloc(256, sizeof(struct MYANCHOR_V *));
	for(ui = 0; ui < 256; ui++)
	{
		mykmer_v2[ui] = (struct MYKMER_V *)malloc(sizeof(struct MYKMER_V));
		mykmer_v2[ui]->size = mykmer_v2[ui]->cap = 0;
		mykmer_v2[ui]->buffer = NULL;
		mykmer_v2[ui]->i = NULL;
		mykmer_v2[ui]->p = NULL;
		myanchor_vs[ui] = (struct MYANCHOR_V *)malloc(sizeof(struct MYANCHOR_V));
		myanchor_vs[ui]->size = myanchor_vs[ui]->cap = 0;
		myanchor_vs[ui]->buffer = NULL;
		myanchor_v2[ui] = (struct MYANCHOR_V *)malloc(sizeof(struct MYANCHOR_V));
		myanchor_v2[ui]->size = myanchor_v2[ui]->cap = 0;
		myanchor_v2[ui]->buffer = NULL;
	}
	//
	clock_t start0, end0, sme0 = 0;
	double usersys0;
	start0 = clock();
	fr = open_all_filereader(1, reffn, 0);
	free(reffn);
	seq = init_biosequence();
	bi = 0;
	while(readseq_filereader(fr, seq))
	{
		mystep1(seq, ksize, wsize, bi, mykmer_vs, nucleobase, mask);
		bi++;
	}
	free_biosequence(seq);
	close_filereader(fr);
	for(ui = 0; ui < mask[5]; ui++)
	{
		mykmer_vs[ui]->i = (u4i *)calloc(mykmer_vs[ui]->size, sizeof(u4i));
		mykmer_vs[ui]->p = (u4i *)calloc(mykmer_vs[ui]->size, sizeof(u4i));
		mystep2(mykmer_vs[ui], mykmer_v2, mask);
	}
	end0 = clock();
	sme0 += end0 - start0;
	usersys0 = (double)sme0 / CLOCKS_PER_SEC;
	fprintf(stderr, "new index = %f sec.\n", usersys0);
	//
	outfp = fopen(outfn, "w");
	//clock_t start, end, sme = 0, start1, end1, sme1 = 0, start2, end2, sme2 = 0, start3, end3, sme3 = 0;
	clock_t start, end, sme = 0;
	//double usersys, usersys1, usersys2, usersys3;
	double usersys;
	start = clock();
	fr = open_all_filereader(1, qryfn, 0);
	free(qryfn);
	seq = init_biosequence();
	//start1 = clock();
	while(readseq_filereader(fr, seq))
	{
		mystep3(seq, ksize, wsize, mykmer_v2[0], nucleobase, mask);
		for(bi = 1; bi < 10; bi++)
		{
			readseq_filereader(fr, seq);
			mystep3(seq, ksize, wsize, mykmer_v2[0], nucleobase, mask);
		}
		//end1 = clock();
		//sme1 += end1 - start1;
		//start2 = clock();
		mystep4(mykmer_v2[0], kmax, myanchor_vs, mykmer_vs, mask);
		//end2 = clock();
		//sme2 += end2 - start2;
		//start3 = clock();
		mystep5(myanchor_vs, bsize, bmin, mmin, outfp, myanchor_v2, seq, mask);
		//end3 = clock();
		//sme3 += end3 - start3;
		//start1 = clock();
	}
	//end1 = clock();
	//sme1 += end1 - start1;
	free_biosequence(seq);
	close_filereader(fr);
	end = clock();
	sme += end - start;
	//usersys1 = (double)sme1 / CLOCKS_PER_SEC;
	//usersys2 = (double)sme2 / CLOCKS_PER_SEC;
	//usersys3 = (double)sme3 / CLOCKS_PER_SEC;
	usersys = (double)sme / CLOCKS_PER_SEC;
	//fprintf(stderr, "step1 = %f sec.\n", usersys1);
	//fprintf(stderr, "step2 = %f sec.\n", usersys2);
	//fprintf(stderr, "step3 = %f sec.\n", usersys3);
	fprintf(stderr, "total = %f sec.\n", usersys);
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
	for(ui = 0; ui < 256; ui++)
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
	free(mykmer_v2);
	free(myanchor_vs);
	free(myanchor_v2);
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
	" cnt      cnt a file\n"
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

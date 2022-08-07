#include "TodMapper.h"

int usage_idx()
{
	fprintf(stdout,
        "usage: TodMapper idx [options] <input file>\n"
        " -o <string> index file, [NULL]\n"
        " -k <int>    kmer size, [15]\n"
        " -w <int>    window size, [10]\n"
        " -v          verbose\n"
        "\n"
        "for example: index human reference genome GRCh38p13.fa\n"
        "TodMapper idx -o GRCh38p13.idx GRCh38p13.fa\n"
        );
        return 1;
}

int main_idx(int argc, char **argv)
{
	FileReader *fr;
        BioSequence *seq;
	FILE *fp;
	char *idxfn;
	b4i c, ksize, wsize, ver, pre_idx, min_idx, bi, bj, bk, bl;
	u1i nucleobase[256];
	u8i key[4], mask[4], ui, uj;
        //u8i key[4], mask[4], ui;
	struct MINIMIZER_T pre_elm[64], min_elm;
	struct MINIMIZER_V *minimizer_v;
	//
	idxfn = NULL;
	ksize = 15;
	wsize = 10;
        ver = 0;
        while((c = getopt(argc, argv, "ho:k:w:v")) != -1)
        {
                switch(c)
                {
			case 'o':
                                idxfn = optarg;
                                break;
			case 'k':
                                ksize = atoi(optarg);
                                break;
			case 'w':
                                wsize = atoi(optarg);
                                break;
                        case 'v':
                                ver++;
                                break;
                        default:
                                return usage_idx();
                }
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
	minimizer_v = (struct MINIMIZER_V *)malloc(sizeof(struct MINIMIZER_V));
	minimizer_v->size = minimizer_v->cap = 0;
	//
	if(optind < argc)
        {
                fr = open_all_filereader(1, argv + argc - 1, 0);
        }
        else
        {
                return usage_idx();
        }
        seq = init_biosequence();
	bi = 0;
	//while(bi < 3 && readseq_filereader(fr, seq))
        while(readseq_filereader(fr, seq))
	{
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
                                        //fprintf(stdout, "%llu %d %d\n", pre_elm[pre_idx].s, bi, bj - ksize + 1);
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
						/*for(bl = 0; bl < wsize; bl++)
						{
							if(pre_elm[bl].s == min_elm.s)
							{
								minimizer_push(minimizer_v, pre_elm[bl]);
							}
						}*/
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
								/*for(bl = 0; bl < wsize; bl++)
								{
									if(pre_elm[bl].s == min_elm.s)
									{
										minimizer_push(minimizer_v, pre_elm[bl]);
									}
								}*/
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
					//fprintf(stdout, "%d %d %d\n", -1, bi, bj - ksize + 1);
				}
			}
			else
			{
				bk = 0;
				min_elm.s = UINT64_MAX;
				//fprintf(stdout, "%d %d %d\n", -1, bi, bj - ksize + 1);
			}
		}
		bi++;
	}
	free_biosequence(seq);
        close_filereader(fr);
	//
	if(idxfn)
	{
		fp = fopen(idxfn, "w");
		/*for(ui = 0; ui < minimizer_v->size; ui++)
		{
			for(bi = 56; bi >= 0; bi = bi - 8)
			{
				nucleobase[128] = (minimizer_v->buffer[ui].s >> bi) & 255;
				fprintf(fp, "%c", nucleobase[128]);
			}
			for(bi = 56; bi >= 0; bi = bi - 8)
                        {
                                nucleobase[128] = (minimizer_v->buffer[ui].ip >> bi) & 255;
                                fprintf(fp, "%c", nucleobase[128]);
                        }
		}*/
		for(ui = 0; ui < minimizer_v->size; ui++)
                {
			uj = minimizer_v->buffer[ui].s;
			fprintf(fp, "%llu ", uj);
			uj = (minimizer_v->buffer[ui].ip >> 32) & mask[3];
			fprintf(fp, "%llu ", uj);
                        uj = minimizer_v->buffer[ui].ip & mask[3];
                        fprintf(fp, "%llu\n", uj);
                }
		fclose(fp);
	}
	else
	{
		return usage_idx();
	}
	//
	free(minimizer_v->buffer);
	free(minimizer_v);
	//
        return 0;
}

int usage_aln()
{
	fprintf(stdout,
	"usage: TodMapper aln [options] <input file>\n"
	" -o <string> index file, [NULL]\n"
	" -k <int>    kmer size, [15]\n"
	" -w <int>    window size, [10]\n"
	" -b <int>    bin size, [256]\n"
	" -v          verbose\n"
	"\n"
	"for example: align multi sequence qry.fa\n"
	"TodMapper aln -o qry.aln qry.fa\n"
	);
        return 1;
}

int main_aln(int argc, char **argv)
{
	FileReader *fr;
	BioSequence *seq;
	FILE *fp;
	char *alnfn;
	b4i c, ksize, wsize, bsize, ver, bi;
	u8i mask[4], qrybin_size, ui, uj;
	struct MINIMIZER_V *minimizer_v;
	struct QRYBIN_V **qrybin_v;
	//
	alnfn = NULL;
	ksize = 15;
	wsize = 10;
	bsize = 256;
	ver = 0;
	while((c = getopt(argc, argv, "ho:k:w:b:v")) != -1)
	{
		switch(c)
		{
			case 'o':
				alnfn = optarg;
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
			case 'v':
				ver++;
				break;
			default:
				return usage_aln();
		}
	}
	//
	mask[3] = (1ULL << 32) - 1;
	minimizer_v = (struct MINIMIZER_V *)malloc(sizeof(struct MINIMIZER_V));
	minimizer_v->size = minimizer_v->cap = 0;
	minimizer_v->buffer = NULL;
	//
	if(optind < argc)
	{
		fr = open_all_filereader(1, argv + argc - 1, 0);
	}
	else
	{
		return usage_aln();
	}
	seq = init_biosequence();
	bi = 0;
	readseq_filereader(fr, seq);
	qrybin_size = seq->seq->size;
	if(qrybin_size % bsize == 0)
	{
		qrybin_size = qrybin_size - (qrybin_size % bsize);
		qrybin_size = qrybin_size / bsize;
	}
	else
	{
		qrybin_size = qrybin_size - (qrybin_size % bsize);
		qrybin_size = qrybin_size / bsize;
		qrybin_size = qrybin_size + 1;
	}
	qrybin_v = (struct QRYBIN_V **)calloc(qrybin_size, sizeof(struct QRYBIN_V *));
	for(ui = 0; ui < qrybin_size; ui++)
	{
		qrybin_v[ui] = (struct QRYBIN_V *)malloc(sizeof(struct QRYBIN_V));
		qrybin_v[ui]->size = qrybin_v[ui]->cap = 0;
		qrybin_v[ui]->buffer = NULL;
	}
	initqrybin(seq, ksize, wsize, bsize, bi, minimizer_v, qrybin_v);
	bi++;
	//
	if(alnfn)
	{
		fp = fopen(alnfn, "w");
		fprintf(fp, "minimizer_v\n");
		for(ui = 0; ui < minimizer_v->size; ui++)
		{
			uj = minimizer_v->buffer[ui].s;
			fprintf(fp, "%llu ", uj);
			uj = (minimizer_v->buffer[ui].ip >> 32) & mask[3];
			fprintf(fp, "%llu ", uj);
			uj = minimizer_v->buffer[ui].ip & mask[3];
			fprintf(fp, "%llu\n", uj);
		}
		for(ui = 0; ui < qrybin_size; ui++)
		{
			fprintf(fp, "qrybin_v[%llu]\n", ui);
			for(uj = 0; uj < qrybin_v[ui]->size; uj++)
			{
				fprintf(fp, "%llu\n", qrybin_v[ui]->buffer[uj].s);
			}
		}
		fclose(fp);
	}
	else
	{
		return usage_aln();
	}
	//
	while(readseq_filereader(fr, seq))
	{
		updateqrybin(seq, ksize, wsize, bsize, bi, minimizer_v, qrybin_v);
		bi++;
	}
	free_biosequence(seq);
	close_filereader(fr);
	//
	if(minimizer_v->size > 0)
	{
		free(minimizer_v->buffer);
	}
	free(minimizer_v);
	for(ui = 0; ui < qrybin_size; ui++)
	{
		if(qrybin_v[ui]->size > 0)
		{
			free(qrybin_v[ui]->buffer);
		}
		free(qrybin_v[ui]);
	}
	free(qrybin_v);
	//
	if(_DEBUG_LOG_)
	{
	}
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
        " idx      index a file\n"
        " aln      align a file\n"
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
	if(strcasecmp("idx", argv[1]) == 0) return main_idx(argc - 1, argv + 1);
	if(strcasecmp("aln", argv[1]) == 0) return main_aln(argc - 1, argv + 1);
	if(strcasecmp("-h", argv[1]) == 0) return usage();
	if(strcasecmp("--help", argv[1]) == 0) return usage();
	fprintf(stderr, " -- unknown command '%s' -- \n", argv[1]);
	return 1;
}

/*
 * main.cpp
 *
 *  Created on: 2021年7月8日
 *      Author: vm
 */
#include "Basic.h"
#include "Graph.h"
#include "Index.h"
#include "BasicMethod.h"
#include "TopoLevel.h"
#include "genIndexMethods.h"
#include "queryMethods.h"
#include "transf.h"
#include "geneQuery.h"

void print_usage() {
	cout << "minBL Reachability Index" << endl;
	cout << "Usage: minBL -g <GRAPHFILE> -M <METHOD> -query <QUERYFILE> -index <INDEXFILE> [-k <SIZE>] [-d <DELTA>]" << endl;
	cout << "Description:\n"
         "-g <GRAPHFILE> is the name of the input graph.\n"
		 "-M <METHOD> is the method of constructing index or handling query or filtering or comparing.\n"
         "-query <QUERYFILE> contains a line of <src> <dest> <reachability> for each query.\n"
		 "-index <INDEXFILE> is the index of minBL.\n"
         "-k <SIZE> sets the value of parameter k, default value is 2.\n"
         "-d <DELTA> sets the value of parameter delta, default value is 2.\n"
	     "-tl 1 is used to compute topological level or handle query by topological level.\n"
	     "-rand 1 is used to compute the index of minBL or handle query by minBL index.\n"
	     "-randk 1 is used to compute the first k index of minBL or handle query by k-minBL index.\n"
	     "-randk 2 is used to compute the first k index of tfd-minBL.\n"
	     "-randk 3 is used to compute the first k index of ctfd-minBL.\n"
	     "-tfd 1 is used to compute the index of tfd-minBL or handle query by tfd-minBL index.\n"
	     "-tfdd 1 is used to compute the index of ctfd-minBL or handle query by ctfd-minBL index.\n" << endl;
}


int main(int argc, char *argv[])
{
	unsigned label_topolevel=0;
	unsigned label_MinV_rand=0;
	unsigned label_MinV_k=0;
	unsigned label_MinV_tfd=0;
	unsigned label_MinV_tfd_delta=0;
	unsigned label_ip = 0;
	unsigned k=2;
	unsigned delta=2;
	unsigned w=0;

	char M='A';
	string h="";
    char * input_graph_filename;
	char * indexFileName=NULL;
	char * queryFileName=NULL;
    char * ipIndexFileName;
	char * f1=NULL;
	char * f2=NULL;
	char * f3=NULL;
//	char * transf_fileName;
	char * gene_fileName;

	for (int i = 0; i < argc; ++i) {
		string arg = argv[i];
		if (arg == "-g") {         // graph file
			input_graph_filename = argv[++i];
		}else if (arg == "-index") {  // index size constraint
            indexFileName = argv[++i];
        }else if (arg == "-ipIndex") {  // delta-minV BL
            ipIndexFileName = argv[++i];
        }else if (arg == "-geneFileName") {  // index size constraint
            gene_fileName = argv[++i];
        }else if (arg == "-query") {  // index size constraint
            queryFileName = argv[++i];
        }else if (arg == "-f1") {  // index size constraint
            f1 = argv[++i];
        }else if (arg == "-f2") {  // index size constraint
            f2 = argv[++i];
        }else if (arg == "-f3") {  // index size constraint
            f3 = argv[++i];
        }else if (arg == "-k") {  // index size constraint
			k = atoi(argv[++i]);
		}else if (arg == "-d") {  // delta-minV BL
            delta = atoi(argv[++i]);
        }else if (arg == "-cho") {  // index size constraint
			w = atoi(argv[++i]);
		}
		/*else if (arg == "-trans") {  // index size constraint
			transf_fileName = argv[++i];
		}*/
		else if (arg == "-rand") {  // delta-minV BL
			label_MinV_rand = atoi(argv[++i]);
		}else if (arg == "-randk") {  // delta-minV BL
			label_MinV_k = atoi(argv[++i]);
		}else if (arg == "-tfd") {  // delta-minV BL
			label_MinV_tfd = atoi(argv[++i]);
		}else if (arg == "-tfdd") {  // delta-minV BL
			label_MinV_tfd_delta = atoi(argv[++i]);
		}else if (arg == "-tl") {  // delta-minV BL
			label_topolevel = atoi(argv[++i]);
		}else if (arg == "-ip") {  // delta-minV BL
            label_ip = atoi(argv[++i]);
        }else if (arg == "-M") {  // delta-minV BL
			M = argv[++i][0];
		}else if (arg == "-h") {  // delta-minV BL
			h = "help";
		}
	}

	if(h == "help") {
		print_usage();
		return 1;
	}
	if(M=='I')
	{
		//1) load graph
		g = new Graph(input_graph_filename);

		//2) constructing index
		struct indexParaSettings p;
		p.delta=delta;
		p.k=k;
		p.g=g;
		p.input_graph_filename=input_graph_filename;
		p.label_MinV_k=label_MinV_k;
		p.label_MinV_rand=label_MinV_rand;
		p.label_MinV_tfd=label_MinV_tfd;
		p.label_MinV_tfd_delta=label_MinV_tfd_delta;
		p.label_topolevel=label_topolevel;
		genIndex(p);
	}
	else if(M=='Q')
	{
		struct queryParaSettings p;
		p.indexFileName=indexFileName;
		p.k=k;
		p.label_MinV_k=label_MinV_k;
		p.label_MinV_rand=label_MinV_rand;
		p.label_MinV_tfd=label_MinV_tfd;
		p.label_MinV_tfd_delta=label_MinV_tfd_delta;
		p.label_topolevel=label_topolevel;
		p.queryFileName=queryFileName;
		p.input_graph_filename = input_graph_filename;
        MinVLabel_RQ(p);
//		MinVLabel_RQ_third(p);
	}
	else if(M=='F')
	{
        struct queryParaSettings p;
        p.indexFileName=indexFileName;
        p.k=k;
        p.label_MinV_k=label_MinV_k;
        p.label_topolevel=label_topolevel;
        p.label_ip = label_ip;
        p.queryFileName=queryFileName;
        p.input_graph_filename = input_graph_filename;
        p.ipIndexFileName = ipIndexFileName;
        Filter_Test(p);
	}
	else if(M=='C')
	{
		isEqual_queryResult(f1, f2);
	}
	else if(M=='G')
	{
		string filename(gene_fileName);
		randQueryPair(input_graph_filename,1000000,filename);
		struct queryParaSettings p;
		p.indexFileName=indexFileName;
		p.k=k;
		p.label_MinV_k=label_MinV_k;
		p.label_MinV_rand=label_MinV_rand;
		p.label_MinV_tfd=label_MinV_tfd;
		p.label_MinV_tfd_delta=label_MinV_tfd_delta;
		p.label_topolevel=label_topolevel;
		p.queryFileName=gene_fileName;
		p.input_graph_filename = input_graph_filename;
		MinVLabel_makeR(p);
	}
    else if(M=='N') {
        struct queryParaSettings p;
        p.indexFileName = indexFileName;
        p.k = k;
        p.label_MinV_k = label_MinV_k;
        p.label_MinV_rand = label_MinV_rand;
        p.label_MinV_tfd = label_MinV_tfd;
        p.label_MinV_tfd_delta = label_MinV_tfd_delta;
        p.label_topolevel = label_topolevel;
        p.input_graph_filename = input_graph_filename;
        p.queryFileName=gene_fileName;
        negativeQueryPair(p);
    }
    else if(M=='P'){
        string fileName(gene_fileName);
        positiveQueryPair(g, 1000000, fileName);
    }
	else if(M=='W')
	{
		string str(input_graph_filename);
		transfer(str,w);
	}
    else if(M=='T')
    {
        compare_filter(f1, f2, f3);
    }
}


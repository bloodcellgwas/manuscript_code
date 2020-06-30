#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <unistd.h>
// William Astle
int main(int argc, char** argv)
{

	std::cout<<"currently assuming no offset"<<std::endl;	
	std::string input_filename;
	std::vector<std::string> input_filenames_vec;
	std::string output_filename;
	std::ifstream test_ifstream;
	std::ofstream test_ofstream;

  	std::string input_weights;
	std::stringstream weight_stream;
	std::string weight_string;
	std::vector<double> weight_vec(3);
	std::vector< std::vector<double> > weight_vec_vec;	
	std::string weight_label_string;
	std::vector<std::string> weight_label_vec;	

	int c;
	bool intercept=false;
	std::cout<< "wja hacked this and added the line below to force intercept inclusion";
	intercept=true;
	

	size_t n_meta_cols=6;	
//	opterr = 0;
	while ((c = getopt (argc, argv, "g:o:w:i")) != -1)
	switch (c)
     	 {
		case 'g':
        		input_filename.assign(optarg);
			input_filenames_vec.push_back(input_filename);
			test_ifstream.open(optarg);
			if(!test_ifstream.good())
				std::cerr<<"Input file read error"<<std::endl;
			test_ifstream.close();
			break;
   	   	case 'o':
      			output_filename.assign(optarg);
			test_ofstream.open(optarg);

			if(!test_ifstream.good())
                                std::cerr<<"Output file write error"<<std::endl;
                        test_ifstream.close();
			break;
     		case 'w':
			input_weights.assign(optarg);
			weight_stream.str(input_weights);
			getline(weight_stream, weight_label_string, ':');
			weight_label_vec.push_back(weight_label_string);
			
			getline(weight_stream, weight_string, ',');
			weight_vec[0]=atof(weight_string.c_str());
        		getline(weight_stream, weight_string, ',');
			weight_vec[1]=atof(weight_string.c_str());
        		getline(weight_stream, weight_string);
			weight_vec[2]=atof(weight_string.c_str());
        		weight_vec_vec.push_back(weight_vec);
			weight_stream.clear();
			break;
		case 'i':
			intercept=true;
		case '?':
			if((optopt!='g')&(optopt!='o')&(optopt!='w')&(optopt!='i'))
				abort();
			break; 
		default:
      			abort ();
      }
 
      if(input_filenames_vec.size()==0)
	{
		std::cerr<<"No input Oxford genotype files specified"<<std::endl;
		abort();
	}
      if(output_filename.length()==0)
	{
		std::cerr<<"No file output file specified"<<std::endl;
		abort();
	}

	if(weight_label_vec.size()==0)
	{
		weight_label_vec.push_back("add");
		weight_vec[0]=-1.0;
		weight_vec[1]=0.0;
		weight_vec[2]=1.0;
		weight_vec_vec.push_back(weight_vec);
		weight_label_vec.push_back("dom");
		weight_vec[0]=0.0;
		weight_vec[1]=1.0;
		weight_vec[2]=0.0;
		weight_vec_vec.push_back(weight_vec);
	
	}
	

	size_t no_datasets=weight_label_vec.size();
	std::ifstream input_file_count_n;
	std::string line;
	std::string chr;
	std::string id;
	std::string rsid;
	size_t bp;
	std::string ref;
	std::string alt;
	double prob, probAA, probRA, probRR;
	
	std::vector<size_t> cum_n_vec(input_filenames_vec.size()+1,0);
	
	size_t file_num=0;


	size_t n,p;

	
	for(std::vector<std::string>::iterator it = input_filenames_vec.begin(); it != input_filenames_vec.end(); ++it)
	{
		size_t local_p=1;
		input_file_count_n.open((*it).c_str());
		if(!input_file_count_n.good())
		{
			std::cerr<<"Input file read error"<<std::endl;
			abort();	
		}
		std::cout<<"Determining n for file "<<*it<<"... ";
		std::getline(input_file_count_n, line);	

		size_t n_spaces=0;	
		for (size_t i = 0; i < line.size(); ++i)
			if (line[i] == ' ') ++n_spaces;
		
		for(size_t i=file_num;i<input_filenames_vec.size();i++)
		{	
			cum_n_vec[i+1]+=(n_spaces-n_meta_cols+1)/3;
		}
		while(std::getline(input_file_count_n, line))
			local_p++;
	
		input_file_count_n.close();

		std::cout<<"n="<<(n_spaces-n_meta_cols+1)/3<<std::endl;
		
		if(file_num==0)
			p=local_p;
		if(p!=local_p)
		{
			std::cerr<<"Input files have different numbers of variants"<<std::endl;
			abort();
		}
		input_file_count_n.close();
		file_num++;
	}
	std::cout<<p<<" variants in all files"<<std::endl;

	size_t total_n=cum_n_vec[cum_n_vec.size()-1];
	std::cout<<"n summed over all files is "<<total_n<<std::endl;
	// allocate mem
	size_t data_size;

	data_size=static_cast<size_t>(total_n*(p+static_cast<size_t>(intercept)));

	std::vector<double* > pdata_vec(no_datasets,NULL);
	for(size_t ds=0;ds<no_datasets;ds++)
	{
		pdata_vec[ds] = new double[data_size];

		if(intercept)
			for(size_t it=0;it<total_n;it++)
				pdata_vec[ds][it]=1.0;
	}
	char delim[]=" ";
	char* token=NULL;
	
	std::ifstream input_file;
	size_t it_p; 
	size_t it_n;
	size_t space_count=0;
	char* unconverted;
	size_t ds;
	double data_value;

	file_num=0;	
	for(std::vector<std::string>::iterator it = input_filenames_vec.begin(); it != input_filenames_vec.end(); ++it)
	{
		input_file.open((*it).c_str());
		if(!input_file.good())
		{
			std::cerr<<"Input file read error"<<std::endl;
			abort();	
		}
		it_p=0;
		while(std::getline(input_file, line))
		{
			char* line_str = new char[line.length()+1];
			std::strcpy (line_str, line.c_str());
			space_count=0;
			it_n=0;
			
			for(token=strtok(line_str, delim); token!=NULL; token=strtok(NULL, delim))
			{
				if(space_count+1>n_meta_cols)
				{
					data_value=strtod(token, &unconverted);

					if(space_count%3==0)
						for(ds=0;ds<no_datasets;ds++)
							pdata_vec[ds][(it_p+static_cast<size_t>(intercept))*total_n+cum_n_vec[file_num]+it_n]=weight_vec_vec[ds][0]*data_value;
	
					if(space_count%3==1)
						for(ds=0;ds<no_datasets;ds++)
							pdata_vec[ds][(it_p+static_cast<size_t>(intercept))*total_n+cum_n_vec[file_num]+it_n]+=weight_vec_vec[ds][1]*data_value;
					if(space_count%3==2)
					{	
						for(ds=0;ds<no_datasets;ds++)
							pdata_vec[ds][(it_p+static_cast<size_t>(intercept))*total_n+cum_n_vec[file_num]+it_n]+=weight_vec_vec[ds][2]*data_value;
						it_n++;
					}	
				}
				space_count++;
			}
			delete[] line_str;
			if(!((it_p+1)%20))	
				std::cout<<"Processing line "<<it_p+1<<" of "<<p<<" in file "<<*it<<"\r";
			it_p++;
		}
		std::cout<<"Processing line "<<it_p<<" of "<<p<<" in file "<<*it<<std::endl;
		input_file.close();
		file_num++;
	}
	

	std::cout<<"Writing to h5..."<<std::endl;
			
	hid_t       file, dataset;         /* file and dataset handles */
	hid_t       dataspace;   /* handles */
	
	hsize_t     dimsf[2];              /* dataset dimensions */
	herr_t      status;                             

	file = H5Fcreate(output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

 	dimsf[0] = p+static_cast<size_t>(intercept);
    	dimsf[1] = total_n;

	for(ds=0;ds<no_datasets;ds++)
	{
		dataspace = H5Screate_simple(2, dimsf, NULL); 
		dataset = H5Dcreate1(file, (weight_label_vec[ds]).c_str(), H5T_IEEE_F32BE, dataspace, H5P_DEFAULT);		
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pdata_vec[ds]);
		H5Sclose(dataspace);	
		H5Dclose(dataset);
	}	
	H5Fclose(file);
	
	for(ds=0;ds<no_datasets;ds++)
		delete[] pdata_vec[ds];	


}

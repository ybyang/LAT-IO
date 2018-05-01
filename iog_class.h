
#ifndef IO_GENERAL_CL_H
#define IO_GENERAL_CL_H

#define _CPP_
#include	"io_general.h"
#include        "lat-io.h"
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

#define INF -65535

inline const char** xqcd_type_str()
{
   static const char	*xqcd_type_dim_desc[N_MAX_DIMENSION_TYPES]={"other", "x", "y", "z", "t", "d", "c", "d2", "c2", "complex", "mass", "smear", "displacement", "s_01", "s_02", "s_03", "s_11", "s_12", "s_13", "d_01", "d_02", "d_03", "d_11", "d_12", "d_13", "conf", "operator", "momentum", "direction", "t2", "mass2", "column", "row", "temporary", "temporary2", "temporary3", "temporary4", "errorbar", "operator2", "param", "fit_left", "fit_right", "jackknife", "jackknife2", "jackknife3", "jackknife4", "summary", "channel", "channel2", "eigen", "d_row", "d_col", "c_row", "c_col", "parity", "noise", "evenodd", "disp_x", "disp_y", "disp_z", "disp_t", "t3", "t4", "t_source", "t_current", "t_sink","bootstrap", "nothing"};
   return xqcd_type_dim_desc;
}

template <typename T>
T*    make_endian_switch_buffer(T *pdata, T n)
{
        T     *buf;
        int     i, j;
        char    *p1, *p2;

        if((buf=(T*)malloc(sizeof(T)*n))==NULL)
        {
                printf("An error occured when allocing memory.\n");
                return NULL;
        }

        for(i=0; i<n; i++)
                for(j=0; j<sizeof(T); j++)
                        ((char*)(buf+i))[j]=((char*)(pdata+i))[sizeof(T)-1-j];

        return buf;
}

struct int_rand
{
    int virtual i_rand(int max) const
    {
      return rand()%max;
    }
};

class datatype
{
public:
    filetype type;
    
    int &ndim;
    char name[500];
     int size;
    one_dim *dim;
    
    int name_c2i(char *name);
    void name_i2c(int ind,char *name);
    
    datatype():ndim(type.head.n_dimensions)
//    {type.head.n_dimensions=0;ndim=type.head.n_dimensions;dim=type.head.dimensions;}
    {type.head.n_dimensions=0;dim=type.head.dimensions;}
    
    void add_dimension(int name,int size);

    void add_dimension(int name,int size,int *list);

    void add_dimension(int name,int size,double *list);
    
    void add_dimension(one_dim &src);

    void insert_dimension(int name,int id);
    
    void clear_ind()
    {ndim=0;}
    
    void print()
    {
        char name[100];
        for(int i=0;i<ndim;i++)
        {
          
           std::string name="dim_"+std::string(xqcd_type_str()[dim[i].type]);
           printf("%-20s:(%8d)",name.c_str(),dim[i].n_indices);
           for(int j=0;j<dim[i].n_indices;j++)
           if(dim[i].indices[i]>10000000)
              printf("%12d",dim[i].indices[j]);
           else
              printf("%8d",dim[i].indices[j]);
           printf("\n");        
        }
        printf("\n");
    }
     
     void load_type();
     
     void save_type();
     
     one_dim &fdim(int name)
     {
        return dim[find_dim(name)];
     }
     int nind(int name)
     {
        int ind=find_dim(name);
        if(ind>=0)return dim[ind].n_indices;
        else return 1;
     }
     int find_index(int ind,int name)
     {
        one_dim& dimX=fdim(name);
        for(int i=0;i<dimX.n_indices;i++)
           if(dimX.indices[i]==ind) return i;
        return -1;
     }
     
     void print_void_ind(int name);
     
     int find_dim(int name);
     
     int remove_ind(int name);
     
     int sort_ind(int ind,int *map,bool sort);
     
//     int combine_all_head(*datatype res,int size0);
     
    int common_ind(datatype &res,int name,one_dim &dest);
};

class general_data_base:public datatype
{

     general_data_base* set_pointer(general_data_base &res,general_data_base &tmp);

public:
     
     void set_size(int ind_f,int &size1,int &size2);

//     filetype type;
     double * data;

////////////////////
// initialize
     general_data_base()
     {
         ndim=0;
         data=NULL;
     }
     general_data_base(const char*namex)
     {
         sprintf(name,"%s",namex);
         ndim=0;
         data=NULL;
     }

     void initialize(int flag=0)
     {
        if(data!=NULL){delete [] data;data=NULL;}
        if(ndim>0)
        {  size=1;
           for(int i=0;i<ndim;i++)
             size*=dim[i].n_indices;
          if(flag==0)
          if(size>0)
          { data=new double[size];memset(data,0,size*sizeof(double));}
          else
          {
              printf("the file size is zero\n");
              print(); 
          }
        }     
     }

     general_data_base(filetype &type0,const char * namex)
     {    sprintf(name,"%s",namex);data=NULL;type=type0;initialize(); }
    
     general_data_base(filetype &type0,const char * namex,double* data0)
     {    sprintf(name,"%s",namex);data=NULL;type=type0;initialize();memcpy(data,data0,size*sizeof(double)); }
    
     ~general_data_base(){if(data!=NULL) {delete [] data;data=NULL;}}
// initialize
////////////////////

///////////////
//io
     void clear(){if(data!=NULL&&size>0)memset(data,0,sizeof(double)*size);}
     
     void free_data(){if(data!=NULL){delete [] data;data=NULL;}}
     
     void load();
     
     void load(int name,std::vector<int> &list);
     
     int save()
     {
             FILE    *fp;
             int     count;
             int     i;
             int     *pbufint;
             double  *pbufdouble;
     
             if((fp=fopen(name, "wb"))==NULL)
             {
                     printf("Failed opening file %s for write.\n", name);
                     return  -1;
             }
     
             if(qutils::is_little_endian())
             {
                     if(fwrite(&type, sizeof(filetype), 1, fp)!=1)
                     {
                             printf("Failed writing the head of file %s.\n", name);
                             return  -1;
                     }
             }
             else
             {
                     printf("Converting to little endian\n");
     
                     pbufint=make_endian_switch_buffer<int>((int*)&type, sizeof(filetype)/sizeof(int));
                     if(fwrite(pbufint, sizeof(filetype), 1, fp)!=1)
                     {
                             printf("Failed writing the head of file %s.\n", name);
                             return  -1;
                     }
                     free(pbufint);
             }
     
             count=1;
             for(i=0; i<type.head.n_dimensions; i++)
                     count*=type.head.dimensions[i].n_indices;
             if(qutils::is_little_endian())
             {
                     if(fwrite(data, sizeof(double), count, fp)!=count)
                     {
                             printf("Failed writing the data of file %s.\n", name);
                             return  -1;
                     }
             }
             else
             {
                     pbufdouble=make_endian_switch_buffer<double>(data, count);
                     if(fwrite(pbufdouble, sizeof(double), count, fp)!=count)
                     {
                             printf("Failed writing the data of file %s.\n", name);
                             return  -1;
                     }
                     free(pbufdouble);
             }
     
             fclose(fp);
     
             return 0;
     }

     
     void copy(general_data_base &dest);
     
     template<typename type>
     void gather_file(const char *name, std::vector<type> &list)
     {
         free_data();
         general_data_base file1(""),file2("");
         for(int k=0;k<list.size();k++)
         {
             sprintf(file1.name,name,list[k]);
             FILE *pfile=fopen(file1.name,"rb");
             if(pfile==NULL)continue;
             fclose(pfile);
             file1.load();
             file2.combine(file1,*this);
             file2.copy(*this);
         }
     }
     template<typename type>
     void gather_file_bin(const char *name, std::vector<type> &list,int bin_size=1)
     {
         free_data();
         general_data_base file1(""),file2(""),file3("");
         int icount=0;
         for(int k=0;k<list.size();k++)
         {
             sprintf(file1.name,name,list[k]);
             FILE *pfile=fopen(file1.name,"rb");
             if(pfile==NULL)continue;
             fclose(pfile);
             file1.load();
//             printf("%4d%4d\n",icount%bin_size,icount);
             if(icount%bin_size==0)file1.copy(file3);
             else 
             for(int j=0;j<file3.size;j++)
                file3.data[j]+=file1.data[j];
             
             if(icount%bin_size==bin_size-1)
             { 
//                printf("combine, %4d\n",icount/bin_size);
                for(int j=0;j<file3.size;j++)
                    file3.data[j]/=bin_size;
                file2.combine(file3,*this);
                file2.copy(*this);
             }
             icount++;
         }
     }
     
//io
///////////////
 
     void print_all(bool slim=false);

     void print_all(double cut);

     double * seek(int n,...);
     
     int new_data(general_data_base &dest,int name,int new_size);
     
     int normal();
  
     int jackknife(general_data_base &res);

     int antijack(general_data_base &res);
     
     int sum(general_data_base &res,int name);
     
     int aver(general_data_base &res,int name,int flag=1);
     
     int aver_cfg(general_data_base &res,int aver_flag=1);

     int make_cfg_bin(general_data_base &res,int bin_size);

     int aver_cfg_bin(general_data_base &res,int bin_size);
     
     int nind_cfg();
     
     int pick(general_data_base &res,int name,int ist,int ied=INF);

     int pick(general_data_base &res,int name,int size,int* list,int flag=1);
     
     int pick(general_data_base &res,int name,std::vector<int> list,int flag=0)
     {
        return pick(res,name,list.size(),list.data(),flag);
     }

     int pick(general_data_base &res,int name,std::vector<double> list)
     {
        int *ilist=new int[list.size()];
        for(int i=0;i<list.size();i++) ilist[i]=list[i]*1000000;
        return pick(res,name,list.size(),ilist);
     }
     int pick(general_data_base &res,one_dim &dimX,int type=-1)
     {
        int r_type=(type==-1)?dimX.type:type;
        return pick(res,r_type,dimX.n_indices,dimX.indices,0);
     }

     int exclude(general_data_base &dest,int name,int index);
     
     int mom_aver(general_data_base &dest,int nmom);
     
     int folding(general_data_base &dest,int *fold,int nt2,int it_st);
     
//     int eff_mass(general_data_base &dest,int fold_type,int nt2,int it_st);
     
     int move_ind(general_data_base &dest,int name,int move_flag);

     int move_ind_head(general_data_base &dest,int name,int move_flag);
     
     int combine(general_data_base &a,general_data_base &b,bool sort=true);
     
     int eff_mass(general_data_base &res,int fold_type,int nt2,int it_st);
 
     int bootstrap_sample(general_data_base &res,const int_rand& r=int_rand());

     int bootstrap(general_data_base &res,int nboot,const int_rand& r=int_rand());
     
};

int check_dim(one_dim &a,one_dim &b);

int check_index(datatype &a, datatype &b);

void inline convert_head(filetype &src,LatInfo &dest)
{
     dest.clear();
     dest.resize(src.head.n_dimensions);
     for(int idim=0;idim<dest.size();idim++)
     {
          dest[idim].name=std::string(xqcd_type_str()[src.head.dimensions[idim].type]);
          dest[idim].size=src.head.dimensions[idim].n_indices;
          dest[idim].indices.resize(dest[idim].size);
          char tmp[20];
          for(int ind=0;ind<dest[idim].size;ind++)
          { 
               sprintf(tmp,"%d",src.head.dimensions[idim].indices[ind]);
               dest[idim].indices[ind]=std::string(tmp);
          }
     }
}

void inline convert_head(LatInfo &src,filetype &dest)
{
      dest.head.n_dimensions=src.size();
      for(int idim=0;idim<src.size();idim++)
      {
            for(int i=0;i<N_MAX_DIMENSION_TYPES;i++)
            if(!strcmp(xqcd_type_str()[i],src[idim].name.c_str())||!strcmp(xqcd_type_str()[i],"nothing")){
                dest.head.dimensions[idim].type=i;
                break;
            }
            dest.head.dimensions[idim].n_indices=src[idim].size;
            long num=0,st=0;
            for(int ind=0;ind<src[idim].size;st=0,ind++)
            if(ind<src[idim].indices.size()&&parse_long(num,st,src[idim].indices[ind]))
               dest.head.dimensions[idim].indices[ind]=num;
            else
               dest.head.dimensions[idim].indices[ind]=-ind-1;
      }
}


void inline convert_from_lat(LatData &src,general_data_base &dest)
{
    convert_head(src.info,dest.type);
    dest.initialize();
    memcpy(dest.data,src.res.data(),dest.size*sizeof(double));
}

void inline convert_to_lat(general_data_base &src,LatData &dest)
{
    convert_head(src.type,dest.info);
    dest.res.resize(src.size);
    memcpy(dest.res.data(),src.data,src.size*sizeof(double));
}

#endif

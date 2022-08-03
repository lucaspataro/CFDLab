#include <ctype.h>
#include <cerrno>
#include <cstdio>
#include <vector>
#include "helper.hpp"
#include <fstream>
#include "regex"
#include "iostream"
#include <boost/algorithm/string.hpp>
#include "enums.hpp"
#include "datastructures.hpp"
#include "grid.hpp"
#ifdef gpp9
// gcc Version >= 9
#include "filesystem"
namespace filesystem = std::filesystem;
#else
// gcc Version < 9
#include "experimental/filesystem"
namespace filesystem = std::experimental::filesystem;
#endif // DEBUG

using namespace std;

/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;

int min_int( const int n1, const int n2 )
{
    if( n1 < n2 ) return n1;
    return n2;
}

/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler( int nLine, const char *szFile, const char *szString )
{
    int err = errno;

    fprintf( ERROUT, "%s:%d Error : %s", szFile, nLine, szString );
    fprintf( ERROUT, "\n" );
    
    /* if an error within the c-library occured, an error code can be   */
    /* found in the global variable err                                 */
    if( err != 0 )
    {
	fprintf( ERROUT, "C-Lib   errno    = %d\n", err);
	fprintf( ERROUT, "C-Lib   strerror = %s\n", strerror( err ) );
    }
    exit(1);
}


/*  for comfort */
#define READ_ERROR(szMessage, szVarName, szFileName, nLine) \
  { char szTmp[80]; \
    if( nLine ) \
	sprintf( szTmp, " %s  File: %s   Variable: %s  Line: %d", szMessage, szFileName, szVarName, nLine ); \
    else \
	sprintf( szTmp, " %s  File: %s   Variable: %s ", szMessage, szFileName, szVarName); \
    ERROR( szTmp ); \
  }
    

/* --------------------------------------------------------------------------*/
/* The function searches the datafile fh for the line defining the variable  */
/* szVarName and returns the respctive string including the value of the     */
/* variable. If there's no appropriate line within the datafile, the program */
/* stops with an error messsage.                                             */
/* ATTENTION: The pointer returned refers to a static variable within the    */
/* function. To maintain the string over several program calls, it has to be */
/* copied!!!                                                                 */
/*                                                                           */
char* find_string( const char* szFileName, const char *szVarName )
{ 
    int nLine = 0;
    size_t i;
    FILE *fh = NULL;
    
    static char szBuffer[MAX_LINE_LENGTH];	/* containes the line read  */
                                               /* from the datafile        */

    char* szLine = szBuffer;
    char* szValue = NULL;
    char* szName = NULL;

    /* open file */
    fh = fopen( szFileName, "rt" );
    if( fh == 0 ) 
	READ_ERROR("Could not open file", szVarName, szFileName, 0);

    /* searching */
    while( ! feof(fh) )
    {
	fgets( szLine, MAX_LINE_LENGTH, fh );
	++nLine;

	/* remove comments */
	for( i = 0; i < std::strlen(szLine); i++)
	    if( szLine[i] == '#' )
	    {
		szLine[i] = '\0'; /* Stringende setzen */
		break;
	    }

	/* remove empty lines */
	while( isspace( (int)*szLine ) && *szLine) ++szLine;
	if( strlen( szLine ) == 0) continue; 

	/* now, the name can be extracted */
	szName = szLine;
	szValue = szLine;
	while( (isalnum( (int)*szValue ) || *szValue == '_') && *szValue) ++szValue;
	
	/* is the value for the respective name missing? */
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	*szValue = 0;		/* complete szName! at the right place */
	++szValue;
        
	/* read next line if the correct name wasn't found */
	if( strcmp( szVarName, szName)) continue;

	/* remove all leading blnkets and tabs from the value string  */
	while( isspace( (int)*szValue) ) ++szValue;
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	fclose(fh);
	return szValue;
    }  
   
    READ_ERROR("variable not found", szVarName, szFileName, nLine);
    
    return NULL;		/* dummy to satisfy the compiler  */
} 

void read_string( const char* szFileName, const char* szVarName, char*   pVariable)
{
    char* szValue = NULL;	/* string containg the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as variable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%s", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName,0);

    printf( "File: %s\t\t%s%s= %s\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              pVariable );
}

void read_int( const char* szFileName, const char* szVarName, int* pVariable)
{
    char* szValue = NULL;	/* string containing the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%d", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName, 0);

    printf( "File: %s\t\t%s%s= %d\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              *pVariable );
}

void read_double( const char* szFileName, const char* szVarName, double* pVariable)
{
    char* szValue = NULL;	/* String mit dem eingelesenen Variablenwert */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%lf", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName, 0);

    printf( "File: %s\t\t%s%s= %f\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              *pVariable );
}

/* ----------------------------------------------------------------------- */
/*                   write matrices to a file                              */
/* ----------------------------------------------------------------------- */

void write_matrix( const char* szFileName,       /* filename */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch,		       /* last row */
		 double xlength,	       /* size of the geometry in */
                                               /* x-direction */
		 double ylength,	       /* size of the geometry in */
                                               /* y-direction  */
		   int fFirst ) 	       /* 0 == append, else overwrite*/
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = new float[nSize];
   int k = 0;

   if( fFirst )				/* first call of the function ? */
   {
       fh = fopen( szFileName, "w");	/* overwrite file/write new file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
	   ERROR( szBuff );
       }
       
/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */
   }
   else
   {
       fh = fopen( szFileName ,"a");	/* append to the file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Outputfile %s cannot be opened", szFileName );
	   ERROR( szBuff );
       }
   } 

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
	   tmp[k++] = (float)m[i][j];

   fwrite( tmp, sizeof(float), nSize, fh);

   if( fclose(fh) )
   {
       char szBuff[80];
       sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
       ERROR( szBuff );
   };

   delete[] tmp;
}


void read_matrix( const char* szFileName,       /* filename */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch		       /* last row */
                  ) 	  
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = new float[nSize];
   int k = 0;

       fh = fopen( szFileName, "r");	/* overwrite file/write new file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", szFileName );
	   ERROR( szBuff );
       }


   fread( tmp, sizeof(float), nSize, fh);

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
	   m[i][j]=tmp[k++];

   if( fclose(fh) )
   {
       char szBuff[80];
       /*orig bug:
       sscanf( szBuff, "Inputfile %s cannot be closed", szFileName );*/
       sprintf( szBuff, "Inputfile %s cannot be closed", szFileName );
       ERROR( szBuff );
   };

   free( tmp );
}


/* ----------------------------------------------------------------------- */
/*                      general matrix functions                           */
/* ----------------------------------------------------------------------- */

/*  allocates storage for a matrix                                         */
/*double **matrix( int nrl, int nrh, int ncl, int nch )
{
   int i;
   int nrow = nrh - nrl + 1;	*//* compute number of lines *//*
   int ncol = nch - ncl + 1;	*//* compute number of columns *//*

   double **pArray  = (double **) malloc((size_t)( nrow * sizeof(double*)) );
   double  *pMatrix = (double *)  malloc((size_t)( nrow * ncol * sizeof( double )));

   if( pArray  == 0)  ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  ERROR("Storage cannot be allocated");

   *//* first entry of the array points to the value corrected by the
      beginning of the column *//*
   pArray[0] = pMatrix - ncl;

   *//* compute the remaining array entries *//*
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   // return the value corrected by the beginning of a line
   return pArray - nrl;
}*/


/* deallocates the storage of a matrix  */
void free_matrix( double **m, int nrl, int nrh, int ncl, int nch )
{
   double **pArray  = m + nrl;
   double  *pMatrix = m[nrl]+ncl;

   delete [] pMatrix;
   delete [] pArray;
}

void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
	   m[i][j] = a;
}


/* allocates storage for a matrix */
int **imatrix( int nrl, int nrh, int ncl, int nch )
{
   int i;

   int nrow = nrh - nrl + 1;	/* compute number of rows */
   int ncol = nch - ncl + 1;	/* compute number of columns */

   int **pArray  = (int **) malloc((size_t)( nrow * sizeof( int* )) );
   int  *pMatrix = (int *)  malloc((size_t)( nrow * ncol * sizeof( int )));


   if( pArray  == 0)  ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the
      beginning of the column */
   pArray[0] = pMatrix - ncl;

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}

/* deallocates the storage of a matrix  */
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch )
{
    int **pArray  = m + nrl;
    int  *pMatrix = m[nrl]+ncl;

   
    delete [] pMatrix;
    delete [] pArray;
}

void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
	   m[i][j] = a;
}


matrix<int> read_pgm(Config& config) {
    FILE *input = NULL;
    char line[1024];
    int levels, geo_x, geo_y;
    int **pic = NULL;

    if ((input = fopen(config.geometry.c_str(), "rb")) == 0) {
        char szBuff[80];
        sprintf(szBuff, "Can not read file %s !!!", config.geometry.c_str());
        ERROR(szBuff);
    }

    /* read first line*/
    if(fgets(line, sizeof line, input) == NULL) {
        ERROR("cannot read PGM.");
    }
    if(strncmp("P2", line, 2) != 0) {
        ERROR("file is not a valid PGM file. If file is a valid PGM file, first two characters of the first line must read P2");
    }

    /* skip the comments */
    int is_comment_line = (*line == '#');
    do{
        if(is_comment_line != 0)
	    printf("Skipping comment line: %s", line);
        if(fgets(line, sizeof line, input) == NULL)
            ERROR("cannot read PGM.");
        is_comment_line = (*line == '#');
    }while (is_comment_line);

    /* read the width and height */
    sscanf(line, "%d %d\n", &geo_x, &geo_y);
    printf("Image size: %d x %d\n", geo_x, geo_y);
    
    /* get number of greyscale levels */
    fscanf(input, "%d", &levels);
    
    /* allocate memory for image */
    pic = imatrix(0, geo_x-1, 0, geo_y-1);
    config.total_n_fluid = 0;
    printf("Image initialised...\n");
    for (int j = geo_y-1; j >= 0; --j) {
        for (int i = 0; i < geo_x; ++i) {
            int value;
            fscanf(input, "%d", &value);

            if (value == EOF) {
                fclose(input);
                ERROR("read of geometry file failed!");
            }
            pic[i][j] = value;

            if(pic[i][j] == (int) cell_type::FLUID) {
                config.total_n_fluid += 1;
            }
        }
    }

    // dont scale the boundaries (-2)
    double x_geo_stepper = (geo_x-2.)/config.imax;
    double y_geo_stepper = (geo_y-2.)/config.jmax;


    int x, y, x_l, x_r, y_b, y_t;
    /* allocate memory for scaled image */
    matrix<int> pic_scal(config.l_imax + 2, vector<int>(config.l_jmax + 2, 0)); // imatrix(0, l_imax+2-1, 0, l_jmax+2-1);

    // scale inner points and count fluid cells
    config.n_fluid = 0;
    for(auto j = 0; j < config.l_jmax+2; j++) {
        for(auto i = 0; i < config.l_imax+2; i++) {
            x = config.il - 1 + (int)((1+(i-1) *x_geo_stepper));
            y = config.jb - 1 + (int)((1+(j-1) *y_geo_stepper));
            if (i == 0){
                x = config.il-1;
            }
            if (j ==0) {
                y = config.jb-1;
            }
            if(pic[x][y] == (int) cell_type::FLUID) {
                config.n_fluid += 1;
            }
            pic_scal[i][j] = pic[x][y];
        }
    }


    /* close file */
    fclose(input);
    free_imatrix(pic,0, geo_x-1, 0, geo_y-1);

    return pic_scal;
}

/* ----------------------------------------------------------------------- */
/*                      folder creation functions                          */
/* ----------------------------------------------------------------------- */

/* checks if a respective folder exists */
void check_dir_exists(std::string outputfolder) {
    // Check if the provided string contains a top level directory
    // i.e. "results/output" -> top level directory is results
    filesystem::path folder(outputfolder);

    if (filesystem::is_directory(folder)) {
        std::cout << "Folder: "<< folder <<  " already exists." << std::endl;
    }
    // if the directory is provided like "results", the top level directory doesn't need to be used
    else {
        std::cout << "Folder: "<< folder <<  " doesn't exists." << std::endl;
        if (filesystem::create_directory(folder))
            std::cout << "Folder: "<< folder <<  " created successfully." << std::endl;
    }
}


/**
 * check if a certain file or directory already exists
 * @param path path to the file
 * @return boolean
 */
bool checkIfPathExists(std::string path) {
    std::ifstream f (path.c_str());
    return f.good();
}

/**
 * makes sure that solution path exists
 * @param path path to the solution
 */
void createSolutionPath(std::string path) {
    bool exists = checkIfPathExists(path);
    if (!exists) {
        std::string subPath;
        std::vector<std::string> results;
        boost::split(results, path, [](char c){return c == '/';});
        int length = results.size();
        for (auto i = 0; i < length - 1; i++) {
            subPath.append(results[i]);
            exists = checkIfPathExists(subPath);
            if(!exists) {
                filesystem::create_directory(subPath);
            }
            subPath.append("/");
        }
    }
}

/* ----------------------------------------------------------------------- */
/*                      string format functions                            */
/* ----------------------------------------------------------------------- */

std::string format_duration(int milliseconds) {
    char buffer[50];
    int minutes = milliseconds / 60 / 1000;
    int seconds = (milliseconds - minutes * 60 * 1000) / 1000;
    int millis = milliseconds - minutes * 60 * 1000 - seconds * 1000;

    sprintf(buffer, "%dm %ds %dms", minutes, seconds, millis);

    return std::string(buffer);
}

/* ----------------------------------------------------------------------- */
/*                      compare file functions                             */
/* ----------------------------------------------------------------------- */
/**
* compare whether two files are the same
* @param p1 file1
* @param p2 file2
* @return boolean
*/
bool compareFiles(const std::string& p1, const std::string& p2) {
    std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
}

bool compareFilesEpsilon(const std::string& p1, const std::string& p2, double epsilon) {
    std::ifstream f1(p1);
    std::ifstream f2(p2);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    std::string line_f1;
    std::string line_f2;
    double value_f1;
    double value_f2;
    int count = 0;
    double res = 0;
    bool isPressure = false;

    //iterate through the files
    while (std::getline(f1,line_f1) && std::getline(f2,line_f2)){

        // used for spliting the strings
        istringstream iss_f1(line_f1);
        istringstream iss_f2(line_f2);
        do
        {
            string substring_f1;
            string substring_f2;
            iss_f1 >> substring_f1;
            iss_f2 >> substring_f2;
            if(substring_f1 == "pressure"){
                isPressure = true;
            }
            if(substring_f1 == "SCALARS"){
                isPressure = false;
            }
            //convert to double if possible, then compare.
            if(isPressure) {
                try
                {
                    value_f1 = std::stod(substring_f1);
                    value_f2 = std::stod(substring_f2);
                    res += (value_f1 - value_f2) * (value_f1 - value_f2);
                    count++;
                }
                catch(std::exception& e)
                {
                    //std::cout << "Could not convert string to double: " << substring_f1 << std::endl;
                }
            }
        } while (iss_f1 && iss_f2);
    }
    res /= (count-1); // -1 because of the 1 in "SCALARS pressure float 1"
    res = sqrt(res);
    return res < epsilon;
}

bool compareFilesMinimum(const std::string& p1, const std::string& p2, double epsilon) {
    std::ifstream f1(p1);
    std::ifstream f2(p2);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    std::string line_f1;
    std::string line_f2;
    double value_f1;
    double value_f2;
    int count = 0;
    double res = 0;
    bool isPressure = false;

    vector<double> pressure_f1;
    vector<double> pressure_f2;

    double min_f1 = 0;
    double min_f2 = 0;

    //iterate through the files
    while (std::getline(f1,line_f1) && std::getline(f2,line_f2)){

        // used for spliting the strings
        istringstream iss_f1(line_f1);
        istringstream iss_f2(line_f2);
        do
        {
            string substring_f1;
            string substring_f2;
            iss_f1 >> substring_f1;
            iss_f2 >> substring_f2;
            if(substring_f1 == "pressure"){
                isPressure = true;
            }
            if(substring_f1 == "SCALARS"){
                isPressure = false;
            }
            //convert to double if possible, then compare.
            if(isPressure) {
                try
                {
                    value_f1 = std::stod(substring_f1);
                    value_f2 = std::stod(substring_f2);

                    if (value_f1 < min_f1){
                        min_f1 = value_f1;
                    }
                    if (value_f2 < min_f2){
                        min_f2 = value_f2;
                    }

                    pressure_f1.push_back(value_f1);
                    pressure_f2.push_back(value_f2);

                }
                catch(std::exception& e)
                {
                    //std::cout << "Could not convert string to double: " << substring_f1 << std::endl;
                }
            }
        } while (iss_f1 && iss_f2);
    }


    for(int i=0; i < pressure_f1.size(); i++){
        value_f1 = pressure_f1[i] - min_f1;
        value_f2 = pressure_f2[i] - min_f2;
        res += (value_f1 - value_f2) * (value_f1 - value_f2);
        count++;
    }
    res /= (count-1); // -1 because of the 1 in "SCALARS pressure float 1"
    res = sqrt(res);
    return res < epsilon;
}

/* ----------------------------------------------------------------------- */
/*                      cmd functions                             */
/* ----------------------------------------------------------------------- */
/**
* check for command line argument and return value
*
* @param argn Number of command line arguments
* @param args All command line arguments
* @param argument The command line argument to look up
* @return boolean
*/
std::string get_cmd_argument(int argn, char** args, const std::string& argument) {
    if (argn == 1) {
        // if only argument return empty string since 
        // the first argument is always the programm name
        return std::string();
    }

    // skip first entry
    for(auto i=1; i < argn; i++) {
        // check if argument is the one we are looking for
        // and if a value is present
        if (args[i] == argument && (i + 1) < argn) {
            return std::string(args[i + 1]);
            // skip next entry
            i++;
        }
    }

    return std::string();
}


/* -------------------------------------------------------------- */
/*                       ouput helper                             */
/* -------------------------------------------------------------- */
// template<class T>
void print_matrix(matrix<cell_type> m, int imax, int jmax) {
    for(auto j = jmax + 1; j > -1; j--) {
        for(auto i = 0; i < imax + 2; ++i) {
            printf("%d ",  (int) m[i][j]);
        }
        printf("\n");
    }
}

void print_matrix(cell_type* m, int imax, int jmax) {
    for(auto j = jmax + 1; j > -1; j--) {
        for(auto i = 0; i < imax + 2; ++i) {
            printf("%d ",  (int) m[j*(jmax+2)+i]);
        }
        printf("\n");
    }
}

void print_matrix(matrix<double> m, int imax, int jmax) {
    for(auto j = jmax + 1; j > -1; j--) {
        for(auto i = 0; i < imax + 2; ++i) {
            printf("%f ",  m[i][j]);
        }
        printf("\n");
    }
}

void print_matrix(matrix<int> m, int imax, int jmax) {
    for(auto j = jmax + 1; j > -1; j--) {
        for(auto i = 0; i < imax + 2; ++i) {
            printf("%d",  m[i][j]);
        }
        printf("\n");
    }
}

void print_matrix(MatrixXd m, int imax, int jmax) {
    for(auto j = jmax + 1; j > -1; j--) {
        for(auto i = 0; i < imax + 2; ++i) {
            printf("%f ",  m(i,j));
        }
        printf("\n");
    }
}
void print_domainMat(MatrixXd m) {
    for(auto j = m.cols() -1; j > -1; j--) {
        for(auto i = 0; i < m.rows(); ++i) {
            printf("%f ",  m(i,j));
        }
        printf("\n");
    }
}

void print_grid_types(Grid& grid, int imax, int jmax) {
    for(auto j = jmax + 1; j > -1; j--) {
        for(auto i = 0; i < imax + 2; ++i) {
            printf("%d ",  (int) grid.cell(i,j).type());
        }
        printf("\n");
    }
}

void writeMatrixToFile(std::string fname, MatrixXd & M) {
    std::ofstream myfile(fname);
    for (int j = 0; j < M.cols(); j++) {
        for (int i= 0; i < M.rows()-1; i++) {
            myfile << M(i,j) << ",";
        }
        myfile << M(M.rows()-1,j) << "\n";
    }
    myfile.close();
}

void print_matrix(MatrixXd m) {
    for(auto i = 0; i < m.rows(); ++i) {
        for(auto j = 0; j < m.cols(); ++j) {
            if(m(i, j) != 0) {
                printf("\033[0;31m");
                printf("%2.4f ", m(i, j));
                printf("\033[0m");
            } else {
                printf("%2.0f ", m(i, j));
            }
        }
        printf("\n");
    }
}

void print_config(Config& config) {
    std::cout << "----------------------" << std::endl;
    std::cout << "rank: " << config.rank << std::endl;
    std::cout << "xlength: " << config.xlength << std::endl;
    std::cout << "ylength: " << config.ylength << std::endl;
    std::cout << "Re: " << config.Re << std::endl;
    std::cout << "itermax: " << config.itermax << std::endl;
    std::cout << "imax: " << config.imax << std::endl;
    std::cout << "jmax: " << config.jmax << std::endl;
    std::cout << "l_imax: " << config.l_imax << std::endl;
    std::cout << "l_jmax: " << config.l_jmax << std::endl;
    std::cout << "calcTemp: " << config.calcTemp << std::endl;
    std::cout << "iproc: " << config.iproc << std::endl;
    std::cout << "jproc: " << config.jproc << std::endl;
    std::cout << "levels: " << config.levels << std::endl;
    std::cout << "n_fluid: " << config.n_fluid << std::endl;
    std::cout << "boundary_size: " << config.boundary_size << std::endl;
    std::cout << "num_proc: " << config.num_proc << std::endl;
    std::cout << "rank_l: " << config.rank_l << std::endl;
    std::cout << "rank_r: " << config.rank_r << std::endl;
    std::cout << "rank_t: " << config.rank_t << std::endl;
    std::cout << "rank_b: " << config.rank_b << std::endl;
    std::cout << "il: " << config.il << std::endl;
    std::cout << "ir: " << config.ir << std::endl;
    std::cout << "jb: " << config.jb << std::endl;
    std::cout << "jt: " << config.jt << std::endl;
    std::cout << "omg_i: " << config.omg_i << std::endl;
    std::cout << "omg_j: " << config.omg_j << std::endl;
    std::cout << "t_end: " << config.t_end << std::endl;
    std::cout << "dt: " << config.dt << std::endl;
    std::cout << "omg: " << config.omg << std::endl;
    std::cout << "eps: " << config.eps << std::endl;
    std::cout << "tau: " << config.tau << std::endl;
    std::cout << "alpha: " << config.alpha << std::endl;
    std::cout << "dt_value: " << config.dt_value << std::endl;
    std::cout << "UI: " << config.UI << std::endl;
    std::cout << "VI: " << config.VI << std::endl;
    std::cout << "GX: " << config.GX << std::endl;
    std::cout << "GY: " << config.GY << std::endl;
    std::cout << "PI: " << config.PI << std::endl;
    std::cout << "PR: " << config.PR << std::endl;
    std::cout << "TI: " << config.TI << std::endl;
    std::cout << "T_h: " << config.T_h << std::endl;
    std::cout << "T_c: " << config.T_c << std::endl;
    std::cout << "beta: " << config.beta << std::endl;
    std::cout << "dx: " << config.dx << std::endl;
    std::cout << "dy: " << config.dy << std::endl;
    std::cout << "----------------------" << std::endl;
}
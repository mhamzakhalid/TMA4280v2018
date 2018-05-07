#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <mpi.h>

using namespace std;


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); 
    int size,rank; 
    MPI_Datatype stype; 
     
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n;

    if(rank == 0){
        n=8;
    }

    MPI_Bcast(&n,1,MPI_INT, 0, MPI_COMM_WORLD);

//**************** Generating the mesh *********************//
    int EPP = n/size; 			  // Elements per processor
    double localtopology[EPP][2]={};  	  //Connectivity Matrix "C(K,V)"   
    double v_local[size][EPP+1]={};       //Vector array with information of ghost entries
    double v_global[size][EPP+1]={};
    double meshstep=1./n;
    double localgeometry[n+1];  	  //Setting Geometry  
    for(int i=rank*EPP;i<rank*EPP+EPP;i++){
      
	localtopology[i][0]=i;
	localtopology[i][1]=i+1; 
        localgeometry[i]=i*meshstep;
      
// Matrix that provides information of shared and ghost
     if(rank==0){
                for(int j=0;j<=EPP;j++){
			v_local[rank][j]=j;             
					}
                    } 
      else{  
             for(int j=0;j<EPP;j++){v_local[rank][j]=EPP*rank+j+1;}
	 		v_local[rank][EPP]=EPP*rank;	
	  }
      fflush(stdout); 
		}   

     MPI_Reduce(&v_local, &v_global, size*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     

//******************************* Printing Connectivity on each processor****************************//
   if(rank==0){

	cout<<endl<<"*** This is how local topology looks like for each process ***"<<endl<<endl;
    for(int i=rank*EPP;i<rank*EPP+EPP;i++){
        cout<<"Element: ("<<i<<")  Vertices:  "<<localtopology[i][0]<<"  "<<localtopology[i][1]<<"  ";
					cout<<endl;   }


//******************************* Printing Vector Class that provides information ********************************//

   //if(rank==0){
                cout <<endl<<"*** Printing vector class, with information of shared entries ***"<<endl<<endl;
		for (int i=0;i< size;i++){
        		for(int j=0;j<=EPP;j++){ cout<<"  "<<v_global[i][j]<<" ";}
			cout<<"   ***Process "<<i<<"   Shared: "<<v_global[i][EPP]<<endl;
		}		     		
     }			

   cout<<endl;
    MPI_Finalize();
}

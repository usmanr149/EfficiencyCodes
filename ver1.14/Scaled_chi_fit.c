

int DescendingOrder(){
  int i,j;
  long double en=0,eff=0,err=0;

  //printf("Count from DescendingOrder = %d\n", count0);

  for(i=0;i<count0;i++)
    for(j=i+1;j<count0;j++)
      if (xi0[i]>xi0[j]){
	  en=xi0[i];
	  eff=yi0[i];
	  err=ei0[i];
	  xi0[i]=xi0[j];
	  yi0[i]=yi0[j];
	  ei0[i]=ei0[j];
	  xi0[j]=en;
	  yi0[j]=eff;
	  ei0[j]=err;
      }
}

int main(int argc, char *argv[]){

  int i,j;
  long double bMatrix[ROWS][COLS];
  long double energy, chi;
  long double qMatrix[ROWS][COLS], qVector[ROWS];

  if(argc!=4)
    {
      printf("\n Scaled_chi_fit Filename_with_all data outfilename_for_CHISquared outfilename_for_scaled_eff  \n");
      //exit(-1);
      return 0;
    }

  file_name=argv[1];
  out_file=argv[2];
  out_file2=argv[3];

  remove(out_file);
  remove(out_file2);


  //Check if file opens
  if((inp=fopen(argv[1],"r"))==NULL)
    {
      printf("Cannot open %s.\n",file_name);
      //exit(-1);
      return 0;
    }

   printf("Parameters read from %s.\n",file_name);

   //Read the data and store it in arrays as xi0, yi0 and ei0
    if(fgets(line,132,inp)!=NULL)
      {
	count0=0;
	count=0;
	while(fscanf(inp,"%Lf %Lf %Lf\n",&xi0[count],&yi0[count],&ei0[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi0[count],yi0[count],ei0[count]);
	    count++;
	    count0++;
	  }
      }
    fclose(inp);
    
    //Call the function Matrix where all matrix formation is done
     Matrix();



/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id$";

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>


struct pdb
{
        char recordtype[5],atnam[4],resnam[4];
        int atno,resno,num;
        float x,y,z,r,s,t,time;
        float m,n;
};
struct pdb a;

int
angle5 (float a1x, float a1y, float a1z, float a2x, float a2y, float a2z,
       float a3x, float a3y, float a3z)
{
  float b1x = 0, b1y = 0, b1z = 0, b2x = 0, b2y = 0, b2z = 0,X = 0, Y = 0, ang = 0, m1x = 0, m1y = 0, m1z = 0, n = 0;
  n = sqrt(pow((a2x-a1x),2) + pow((a2y-a1y),2) + pow((a2z-a1z),2));
  b1x = (a1x - a2x)/n;
  b1y = (a1y - a2y)/n;
  b1z = (a1z - a2z)/n;
  n = sqrt(pow((a2x-a3x),2) + pow((a2y-a3y),2) + pow((a2z-a3z),2));
  b2x = (a3x - a2x)/n;
  b2y = (a3y - a2y)/n;
  b2z = (a3z - a2z)/n;
  ang = acos(b1x*b2x + b1y*b2y + b1z*b2z)*180*7/22;
  if (ang >= 70 && ang <= 150)
  {
  return 1;
  }
   else
  {
   return 0;
  }
}
int
angle6 (float a1x, float a1y, float a1z, float a2x, float a2y, float a2z,
       float a3x, float a3y, float a3z)
{
  float b1x = 0, b1y = 0, b1z = 0, b2x = 0, b2y = 0, b2z = 0,X = 0, Y = 0, ang = 0, m1x = 0, m1y = 0, m1z = 0, n = 0;
  n = sqrt(pow((a2x-a1x),2) + pow((a2y-a1y),2) + pow((a2z-a1z),2));
  b1x = (a1x - a2x)/n;
  b1y = (a1y - a2y)/n;
  b1z = (a1z - a2z)/n;
  n = sqrt(pow((a2x-a3x),2) + pow((a2y-a3y),2) + pow((a2z-a3z),2));
  b2x = (a3x - a2x)/n;
  b2y = (a3y - a2y)/n;
  b2z = (a3z - a2z)/n;
  ang = acos(b1x*b2x + b1y*b2y + b1z*b2z)*180*7/22;
  if (ang >= 70 && ang <= 150)
  {
  return 1;
  }
  else
  {
  return 0;
  }
}

int main(int argc,char *argv[])
{
  FILE *fp,*fq,*fd,*fa,*fm,*fb,*fc,*fg;
        int i=0,t=0,atno[17000],resno[17000],num[200],b=6400,framenum=0,residue=0,h,sd;
        int ts,ma,bin=0,o=0,count,count1,count2,count3,count4,count6,count5,ine,d,count10,count11,angfinal,half6cages,final1,che,angfinal1,hj,uc,vc,aj,ac;
      int count7,count8,check,fivememring,check6,sixmemring,half5cages,final,ou[17000],chec,newhalf5cages[17000000];
      float CMg[170000],CMh[170000],CMi[170000];
        float METx, METy, METz, CMx[170000], CMy[170000], CMz[170000],CMd[170000], CMe[170000],CMf[170000],dista;
	  float MTLx, MTLy, MTLz, CMa[170000], CMb[170000], CMc[170000],sumx,sumy,sumz;
	    float x[50000],y[50000],z[50000],distMET,r[1000],s[1000],distMTL,Ox[17000],Oy[17000],Oz[17000],fx[1700000],fy[1700000],fz[1700000];
	      char atnam[1000][4],resnam[1000][4],recordtype[1000][5];
	      float ph,phi_in,ang_in,add,value,qin,ord_in,ord_in1,qin1,ph1,ang_in1,value1,add1;
	      int g,ine1,j,ine2,k,l,ine3,m,c1,ine4,p,q,c2,c3,u[17000],c,e,f,w,v,count9,ine5;
	      float CMj[170000],CMk[170000],CMl[170000],CMm[170000],CMn[170000],CMo[170000],distb,distc,distd,diste,CMp[170000],CMq[170000],CMr[170000];
	      float Ofx1,Ofy1,Ofz1,Ofx2,Ofy2,Ofz2,Ofx3,Ofy3,Ofz3,Ofx4,Ofy4,Ofz4,Ofx5,Ofy5,Ofz5,Ofx6,Ofy6,Ofz6,Ofx7,Ofy7,Ofz7,Ofx8,Ofy8,Ofz8,Ofx9,Ofy9,Ofz9,Ofx10,Ofy10,Ofz10;
	      float Osx1,Osy1,Osz1,Osx2,Osy2,Osz2,Osx3,Osy3,Osz3,Osx4,Osy4,Osz4,Osx5,Osy5,Osz5,Osx6,Osy6,Osz6,Osx7,Osy7,Osz7,Osx8,Osy8,Osz8,Osx9,Osy9,Osz9,Osx10,Osy10,Osz10,Osx11,Osy11,Osz11,Osx12,Osy12,Osz12;
	      float angfive1,ang1,angfive2,ang2,angfive3,ang3,angfive4,ang4,angfive5,ang5,ang6,ang7,ang8,ang9,ang10,ang11;	 
	      float angsix1,angsix2,angsix3,angsix4,angsix5,angsix6;     
  static char *desc[] = {
    "this is a small test program meant to serve as a template ",    "when writing your own analysis tools. The advantage of ",    "using gromacs for this is that you have access to all ",    "information in the topology, and your program will be ",    "able to handle all types of coordinates and trajectory ",    "files supported by gromacs. Go ahead and try it! ",    "This test version just writes the coordinates of an ",    "arbitrary atom to standard out for each frame. You can ",    "select which atom you want to examine with the -n argument."
  };

  static int n=1;

  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    }
  };
  
  t_topology top;
  int        ePBC,ra;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status;
  int        flags = TRX_READ_X;

  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* We don't need any topology information to write the coordinates,
   * but to show how it works we start by writing the name and
   * charge of the selected atom. It returns a boolean telling us
   * whether the topology was found and could be read
   */
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  sfree(xtop);

  n=n-1; /* Our enumeration started on 1, but C starts from 0 */
  /* check that this atom exists */
  if(n<0 || n>(top.atoms.nr)) 
  {
    printf("Error: Atom number %d is out of range.\n",n);
    exit(1);
  }
  
/*	fp=fopen("initial.pdb","rb");
        if(fp==NULL)
                {
                        puts("Cannot Open Source File \n");
                        exit(1);
                }
	while(fp!=NULL)
        {
        fscanf(fp,"%s%d%s%s%d%f%f%f%f%f",&a.recordtype,&a.atno,&a.atnam,&a.resnam,&a.resno,&a.x,&a.y,&a.z,&a.r,&a.s);

        if (feof(fp))
                {
                    break;
                }
	if(strncmp(a.recordtype,"ATOM",4)==0)
		{
		x[i]=a.x; y[i]=a.y; z[i]=a.z; r[i]=a.r; s[i]=a.s;
                atno[i]=a.atno; resno[i]=a.resno;
		for (t=0; t<4; t++)
		recordtype[i][t]=a.recordtype[i][t];
		recordtype[i][4]='\0';
          
		for (t=0; t<3; t++)
			atnam[i][t]=a.atnam[t];
		atnam[i][3]='\0';

		for (t=0; t<3; t++)
			resnam[i][t]=a.resnam[t];
		resnam[i][3]='\0';
		i++;
		}	
	}*/
// fp=fopen("defectcount270-6-500.xvg","w");
  fd=fopen("Noofcages260-6-500run1.xvg","w");
//   fc=fopen("indexdef.xvg","w");
//    fm=fopen("indexOW","w");
//    fa=fopen("indexGu","w");
//    fg=fopen("indexHW","w");
//    fb=fopen("smallcage.pdb","w");
//    fq=fopen("largecage.pdb","w");
//    fa=fopen("../sortfile/cia.txt","r");
//  printf("Atom name: %s\n",*(top.atoms.atomname[n]));
//  printf("Atom charge: %f\n",top.atoms.atom[n].q);
  
  /* The first time we read data is a little special */
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
   printf("i=%d \n",i);
  /* This is the main loop over frames */

            
  do {
    /* coordinates are available in the vector fr.x
     * you can find this and all other structures in
     * the types directory under the gromacs include dir.
     * Note how flags determines wheter to read x/v/f!*/
      
 int s1=-1,s2=-1;	
		    count=0;distMET=0;ine=0;d=0;dista=0;count11=0;half5cages=0;count10=0;u[17000]=0;final=0;ou[17000]=0;chec=0;angfinal=0;angfinal1=0,che=0;final1=0;half6cages=0;vc=0;
count1=0;fx[170000]=0;fy[170000]=0;fz[170000]=0;count4=0;CMx[170000]=0;CMy[17000]=0;CMz[170000]=0;CMd[170000]=0;CMe[170000]=0;CMf[170000]=0;
count6=0;count7=0;count8=0;e=0;distb=0;CMg[170000]=0;CMh[170000]=0;CMi[170000]=0;uc=0;newhalf5cages[17000000]=0;

vc=0;aj=0;ac=0;
angfive1=0;ang1=0;angfive2=0;ang2=0;angfive3=0;ang3=0;angfive4=0;ang4=0;angfive5=0;ang5=0;ang6=0;ang7=0;ang8=0;ang9=0;ang10=0;ang11=0;e=0;sd=0;j=0;hj=0;
Ofx1=0;Ofy1=0;Ofz1=0;Ofx2=0;Ofy2=0;Ofz2=0;Ofx3=0;Ofy3=0;Ofz3=0;Ofx4=0;Ofy4=0;Ofz4=0;Ofx5=0;Ofy5=0;Ofz5=0;Ofx6=0;Ofy6=0;Ofz6=0;Ofx7=0;Ofy7=0;Ofz7=0;Ofx8=0;Ofy8=0;Ofz8=0;Ofx9=0;Ofy9=0;Ofz9=0;Ofx10=0;Ofy10=0;Ofz10=0;
fivememring=0;
angsix1=0;angsix2=0;angsix3=0;angsix4=0;angsix5=0;angsix6=0;
Osx1=0;Osy1=0;Osz1=0;Osx2=0;Osy2=0;Osz2=0;Osx3=0;Osy3=0;Osz3=0;Osx4=0;Osy4=0;Osz4=0;Osx5=0;Osy5=0;Osz5=0;Osx6=0;Osy6=0;Osz6=0;Osx7=0;Osy7=0;Osz7=0;Osx8=0;Osy8=0;Osz8=0;Osx9=0;Osy9=0;Osz9=0;Osx10=0;Osy10=0;Osz10=0;Osx11=0;Osy11=0;Osz11=0;Osx12=0;Osy12=0;Osz12=0;
count2=0;count3=0;count5=0;h=0;ine1=0;ine2=0;k=0;l=0;ine3=0;ine4=0;p=0;q=0;c1=0;c2=0;c3=0;f=0;m=0;h=0;count9=0;diste=0;ine5=0;CMp[17000]=0;CMq[17000]=0;CMr[17000]=0;count10=0;
CMj[170000]=0;CMk[170000]=0;CMl[170000]=0;CMm[170000]=0;CMn[170000]=0;CMo[170000]=0;distb=0;distc=0;distd=0;w=0;v=0;

void allocate_mem(int*** arr, int n, int m)
{
  *arr = (int**)malloc(n*sizeof(int*));
  for(int i=0; i<n; i++)
    (*arr)[i] = (int*)malloc(m*sizeof(int));
     (*arr)[i]=0;
}

//fprintf (fa, "[Guest] \n");
//fprintf (fm, "[OW] \n");
//fprintf (fg, "[HW] \n");


        float distf=0,distg=0;
       int count2=0,count3=0,w1=0,in=0,coun=0,ta[100000][1000]={{0}},ta1[100000][1000]={{0}},ta2[100000][1000]={{0}},ta3[100000][1000]={{0}},gu[10000][1000]={{0}},gu1[10000]={0},gu2[10000]={0},gu3[10000]={0},gu4[10000][1000]={{0}},gu5[10000]={0},gu6[10000]={0},gu7[10000]={0},counts=0,inf=0,inf1=0,counts1=0,gu8[10000]={0},w21=0;
       float distMET=0;
       int index=0,index1=0,index2=0,index3=0,ind,inf2=0,f1=0,u1=0,w6=0,w7=0,w8=0,w9=0,w10=0,w11=0,w12=0,w13=0,w14=0,w15=0,w16=0,w17=0,w18=0,w19=0,w20=0,n1,n2;
      int A1[7][100000] = { {0}, {0} },A2[8][100000] = { {0}, {0} },dr[1000000]={0},dr1[10000]= {0};
     int mem[7][100000]={{0}, {0}},mem1[7][100000]={{0}, {0}};
     int me[8][100000]={{0}},u4=0, u5=0,me1[8][100000]={{0}},trans[100000]={0},trans1[100000]={0};
     int A3[7][100000]= {{0}},A4[100000]= {0},A5[8][100000]= {{0}},A6[8][100000]= {{0}},A7[100000]= {0};


int h = 0, e1 = 0, g = 0,t1=0,t2=0,t3,t4=0,t6=0,max[100]={0},min[100]={0},largecages=0,smallcages=0;
int c=0, f=0,v=0,h1=0,t7=0;
int num[100]={0},in2[100000]={0};
float tim[100]={0};

/*      for (count=0; count<13;count++)
          {
         fscanf(fa,"%f %d",&tim[count],&num[count]);
          }*/

         for(count=0;count<top.atoms.nr;count++)
             {
       //      if (tim[count]==fr.time)
                  {
        //  if (fr.x[count][ZZ] >= 0 && fr.x[count][ZZ] <= 6.16 || fr.x[count][ZZ] >= 13.53 && fr.x[count][ZZ] <= 19.47 )
                   {
            s1 = strcmp (*top.atoms.atomname[count], "Ci");
              s2 = strcmp (*top.atoms.atomname[count], "Ca");
               if (s1==0 || s2==0)
                  {
              w1++;
           in=0; 
                                for (count1=0; count1<top.atoms.nr; count1++)
                                                {
                                        s2 = strcmp (*top.atoms.atomname[count1], "OW1");
                                        s1 = strcmp (*top.atoms.atomname[count1], "OW");
                                           if (s2==0 || s1==0)
                                                  {
  //   distMET = sqrt(pow(((x[count]/10) - fr.x[count1][XX]),2) + pow(((y[count]/10) - fr.x[count1][YY]),2) + pow(((z[count]/10) - fr.x[count1][ZZ]),2));
//  distMET = sqrt(pow((fr.x[num[count]][XX] - fr.x[count1][XX]),2) + pow((fr.x[num[count]][YY] - fr.x[count1][YY]),2) + pow((fr.x[num[count]][ZZ] - fr.x[count1][ZZ]),2));
    distMET = sqrt(pow((fr.x[count][XX] - fr.x[count1][XX]),2) + pow((fr.x[count][YY] - fr.x[count1][YY]),2) + pow((fr.x[count][ZZ] - fr.x[count1][ZZ]),2));
                                        if (distMET<=0.5368 )
                                                {
                                                        in++;
                                               ta[w1][in]=count1;
                                               gu[w1][in]=count;                                                                           
if (in==24)
{
inf1++;
for (coun=1; coun<=in; coun++)
                      {
                      ta1[inf1][coun]=ta[w1][coun];    
                      gu1[inf1]=gu[w1][coun];    
                      in2[inf1]=in;                                        
                      }
}
                                                }
                                                  }
                                                }
                  }
                    }
              }
                  }
        for (coun=1; coun<=inf1; coun++)
                    { 
                      t1=0;t2=0;largecages=0;t3=0,n1=0; 
         for (count2=1; count2<=in2[coun]; count2++)
                      {
                             for (count3=1; count3<=in2[coun]; count3++)
                                                {
     distMET = sqrt(pow((fr.x[ta1[coun][count2]][XX]- fr.x[ta1[coun][count3]][XX]),2) + pow((fr.x[ta1[coun][count2]][YY] -fr.x[ta1[coun][count3]][YY]),2) + pow((fr.x[ta1[coun][count2]][ZZ] - fr.x[ta1[coun][count3]][ZZ]),2));
                                        if (distMET<=0.35 && distMET>0.25)
                                                {
                                                for (count5=1; count5<=in2[coun]; count5++)
                                                            {
         dista = sqrt(pow((fr.x[ta1[coun][count3]][XX] - fr.x[ta1[coun][count5]][XX]),2) + pow((fr.x[ta1[coun][count3]][YY]- fr.x[ta1[coun][count5]][YY]),2) + pow((fr.x[ta1[coun][count3]][ZZ] - fr.x[ta1[coun][count5]][ZZ]),2));
                                    if (dista<= 0.35 && dista>0.25) 
                                    
                                             if (ta1[coun][count2]!=ta1[coun][count5])
                                            {
                                                        for (count6=1; count6<=in2[coun]; count6++)
                                                            {
       distb = sqrt(pow((fr.x[ta1[coun][count5]][XX] - fr.x[ta1[coun][count6]][XX]),2) + pow((fr.x[ta1[coun][count5]][YY]- fr.x[ta1[coun][count6]][YY]),2) + pow((fr.x[ta1[coun][count5]][ZZ] - fr.x[ta1[coun][count6]][ZZ]),2));
                             if (distb<= 0.35 && distb>0.25 )
                             {
                                     if (ta1[coun][count3]!=ta1[coun][count6])
                                     {
                                                        for (count7=1; count7<=in2[coun]; count7++)
                                                            {
      distc = sqrt(pow((fr.x[ta1[coun][count6]][XX] - fr.x[ta1[coun][count7]][XX]),2) + pow((fr.x[ta1[coun][count6]][YY]- fr.x[ta1[coun][count7]][YY]),2) + pow((fr.x[ta1[coun][count6]][ZZ] - fr.x[ta1[coun][count7]][ZZ]),2));
                                  if (distc<= 0.35 && distc>0.25)
                                                   {
                                                          if (ta1[coun][count5]!=ta1[coun][count7])
                                                          {

                                                   for (count8=1; count8<=in2[coun]; count8++)
                                                            {
     distd = sqrt(pow((fr.x[ta1[coun][count7]][XX] - fr.x[ta1[coun][count8]][XX]),2) + pow((fr.x[ta1[coun][count7]][YY]- fr.x[ta1[coun][count8]][YY]),2) + pow((fr.x[ta1[coun][count7]][ZZ] - fr.x[ta1[coun][count8]][ZZ]),2));
                                     if (distd<= 0.35 && distd>0.25 )
                                                   {
                                                           if (ta1[coun][count6]!=ta1[coun][count8])
                                                           {
                                                        if (ta1[coun][count8]==ta1[coun][count2])
                                                        {
                                                           t1++; 
                                                                   inf++;
mem[0][t1] = ta1[coun][count2]; mem[1][t1] =ta1[coun][count3] ; mem[2][t1] = ta1[coun][count5]; mem[3][t1] =ta1[coun][count6]; mem[4][t1] =ta1[coun][count7];
                                                        }
                                                                  for (count9=1; count9<=in2[coun]; count9++)
                                                                {
                       diste = sqrt(pow((fr.x[ta1[coun][count8]][XX] - fr.x[ta1[coun][count9]][XX]),2) + pow((fr.x[ta1[coun][count8]][YY]- fr.x[ta1[coun][count9]][YY]),2) + pow((fr.x[ta1[coun][count8]][ZZ] - fr.x[ta1[coun][count9]][ZZ]),2));
  
                                                                 if (diste<= 0.35 && diste>0.25  )
                                                                 {
                                                                      if (ta1[coun][count9]!=ta1[coun][count7])
                                                                      {

                                                                         if (ta1[coun][count9]==ta1[coun][count2])
                                                                         {
                                                                                 t2++;
                                                                                 inf2++;
                                                                           
                                        me[0][t2] = ta1[coun][count2]; me[1][t2] =ta1[coun][count3] ; me[2][t2] = ta1[coun][count5]; me[3][t2] =ta1[coun][count6]; me[4][t2] =ta1[coun][count7]; me[5][t2] =ta1[coun][count8];

                                                                         }
                                                                     }
                                                                  }
                                                                } 
                                                        
                                                              }
                                                        }
                                                           }

                                                          }
                                                   }
                                                           }
                                      }
                             }
                                                           }
 
                                            }
                                    
                                                          }
                                                 }
                                             
                                                }
                                                }
                                                                    for(h=1;h<=t1;h++)
                                                                           {
                                                                                 check=0;
                                                                         for( t4=h+1;t4<=t1;t4++)
                                                                                  {
                                                                     counts=0;
                                                        for (e1 = 0; e1 <= 4; e1++)
                                                                {
                                                                         for (g = 0; g <=4; g++)
                                                                                  {
                                                                                if (mem[e1][t4]==mem[g][h])
                                                                                                {
                                                                                                    counts++;
                                                                                                }
                                                                                  }
                                                               if (counts >= 5)
                                                                      {
                                                                       check=-1;
                                                                      }
                                                                 }
                                                                       
                                                                                   }
                                                     int check_for6 = 0;
                                                         for (aj = 1; aj <= t2; aj++)
                                                          {
                                                           count11 = 0;
                                                           for (ac = 0; ac <= 5; ac++)
                                                            {
                                                             for (vc = 0; vc <= 4; vc++)
                                                             {
                                                              if (me[ac][aj] == mem[vc][h])
                                                               {
                                                                 count11++;
                                                               }
                                                             }
                                                             if (count11 >= 5)
                                                             {
                                                              check_for6 = -1;
                                                             }
                                                            }
                                                           }
                                                                   if(check == 0 && check_for6==0  )
                                                                       {
                                                       int check1 =  angle5 (fr.x[mem[0][h]][XX], fr.x[mem[0][h]][YY], fr.x[mem[0][h]][ZZ], fr.x[mem[1][h]][XX], fr.x[mem[1][h]][YY], fr.x[mem[1][h]][ZZ],
                                                              fr.x[mem[2][h]][XX], fr.x[mem[2][h]][YY], fr.x[mem[2][h]][ZZ]);
                                                       int check2 =  angle5 (fr.x[mem[1][h]][XX], fr.x[mem[1][h]][YY], fr.x[mem[1][h]][ZZ], fr.x[mem[2][h]][XX], fr.x[mem[2][h]][YY], fr.x[mem[2][h]][ZZ],
                                                                                                                     fr.x[mem[3][h]][XX], fr.x[mem[3][h]][YY], fr.x[mem[3][h]][ZZ]);
                                                       int check3 =  angle5 (fr.x[mem[2][h]][XX], fr.x[mem[2][h]][YY], fr.x[mem[2][h]][ZZ], fr.x[mem[3][h]][XX], fr.x[mem[3][h]][YY], fr.x[mem[3][h]][ZZ],
                                                                            fr.x[mem[4][h]][XX], fr.x[mem[4][h]][YY], fr.x[mem[4][h]][ZZ]);
                                                       int check4 =  angle5 (fr.x[mem[3][h]][XX], fr.x[mem[3][h]][YY], fr.x[mem[3][h]][ZZ], fr.x[mem[4][h]][XX], fr.x[mem[4][h]][YY], fr.x[mem[4][h]][ZZ],
                                                                        fr.x[mem[0][h]][XX], fr.x[mem[0][h]][YY], fr.x[mem[0][h]][ZZ]);
                                                       int check5 =  angle5 (fr.x[mem[4][h]][XX], fr.x[mem[4][h]][YY], fr.x[mem[4][h]][ZZ], fr.x[mem[0][h]][XX], fr.x[mem[0][h]][YY], fr.x[mem[0][h]][ZZ],
                                                                                                                               fr.x[mem[1][h]][XX], fr.x[mem[1][h]][YY], fr.x[mem[1][h]][ZZ]);
                                              //    if (check1 == 1 && check2 == 1 && check3 == 1 && check4 == 1 && check5 == 1 )

                                                                  {
                                                       float d1=sqrt(pow((fr.x[mem[1][h]][XX] - fr.x[mem[4][h]][XX]),2) + pow((fr.x[mem[1][h]][YY]- fr.x[mem[4][h]][YY]),2) + pow((fr.x[mem[1][h]][ZZ] - fr.x[mem[4][h]][ZZ]),2));
                                                       float d4=sqrt(pow((fr.x[mem[2][h]][XX] - fr.x[mem[4][h]][XX]),2) + pow((fr.x[mem[2][h]][YY]- fr.x[mem[4][h]][YY]),2) + pow((fr.x[mem[2][h]][ZZ] - fr.x[mem[4][h]][ZZ]),2));
                                                       float d5=sqrt(pow((fr.x[mem[1][h]][XX] - fr.x[mem[3][h]][XX]),2) + pow((fr.x[mem[1][h]][YY]- fr.x[mem[3][h]][YY]),2) + pow((fr.x[mem[1][h]][ZZ] - fr.x[mem[3][h]][ZZ]),2));
                                                       float d6=sqrt(pow((fr.x[mem[0][h]][XX] - fr.x[mem[3][h]][XX]),2) + pow((fr.x[mem[0][h]][YY]- fr.x[mem[3][h]][YY]),2) + pow((fr.x[mem[0][h]][ZZ] - fr.x[mem[3][h]][ZZ]),2));  
                                                       float d7=sqrt(pow((fr.x[mem[0][h]][XX] - fr.x[mem[2][h]][XX]),2) + pow((fr.x[mem[0][h]][YY]- fr.x[mem[2][h]][YY]),2) + pow((fr.x[mem[0][h]][ZZ] - fr.x[mem[2][h]][ZZ]),2));
                                                           //        if (d1<0.30 || d4<0.30 || d5<0.30 || d6<0.30 || d7<0.30)      
                                                                           {          
                                                                       largecages++;
                                                                   //       fprintf(fd," %d %d %d %d %d %d %d %d \n",mem[0][h],mem[1][h],mem[2][h],mem[3][h],mem[4][h],largecages,coun,gu1[coun]);
                                                                                     
                                                       A1[0][fivememring] = mem[0][h];
                                                       A1[1][fivememring] = mem[1][h];
                                                       A1[2][fivememring] = mem[2][h];
                                                       A1[3][fivememring] = mem[3][h];
                                                       A1[4][fivememring] = mem[4][h];
                                                       A3[0][largecages]=mem[0][h];
                                                       A3[1][largecages]=mem[1][h];
                                                       A3[2][largecages]=mem[2][h];
                                                       A3[3][largecages]=mem[3][h];
                                                       A3[4][largecages]=mem[4][h];
                                                       gu2[largecages]=gu1[coun];
                                                       fivememring++;
                                                     //  index++;
                                                     //  index1++;
                                                  /*   fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[mem[0][h]][XX]*10,fr.x[mem[0][h]][YY]*10,fr.x[mem[0][h]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[mem[1][h]][XX]*10,fr.x[mem[1][h]][YY]*10,fr.x[mem[1][h]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[mem[2][h]][XX]*10,fr.x[mem[2][h]][YY]*10,fr.x[mem[2][h]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[mem[3][h]][XX]*10,fr.x[mem[3][h]][YY]*10,fr.x[mem[3][h]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[mem[4][h]][XX]*10,fr.x[mem[4][h]][YY]*10,fr.x[mem[4][h]][ZZ]*10);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-4,index);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-3,index-4);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-2,index-3);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-1,index-2);
                                                      fprintf(fq,"CONECT%5d%5d\n",index,index-1);  */
                                                                                    
                                                                           }
                                                                  }
                                                              }
                        
                                                                 }
                                                                    for(h1=1;h1<=t2;h1++)
                                                                           {
                                                                                 check6=0;
                                                                         for( c=h1+1;c<=t2;c++)
                                                                                  {
                                                                     counts1=0;
                                                        for (f = 0; f <= 5; f++)
                                                                {
                                                                         for (v = 0; v <=5; v++)
                                                                                  {
                                                                                if (me[f][c]==me[v][h1])
                                                                                                {
                                                                                                    counts1++;
                                                                                                }
                                                                                  }
                                                               if (counts1 >= 6)
                                                                      {
                                                                       check6=-1;
                                                                      }
                                                                 }

                                                                                   }

                                                                             if(check6 == 0  )
                                                                                  {
                                                       int check7 =  angle6 (fr.x[me[0][h1]][XX], fr.x[me[0][h1]][YY], fr.x[me[0][h1]][ZZ], fr.x[me[1][h1]][XX], fr.x[me[1][h1]][YY], fr.x[mem[1][h1]][ZZ],
                                                              fr.x[me[2][h1]][XX], fr.x[me[2][h1]][YY], fr.x[me[2][h1]][ZZ]);
                                                       int check8 =  angle6 (fr.x[me[1][h1]][XX], fr.x[me[1][h1]][YY], fr.x[me[1][h1]][ZZ], fr.x[me[2][h1]][XX], fr.x[me[2][h1]][YY], fr.x[mem[2][h1]][ZZ],
                                                                                                                     fr.x[me[3][h1]][XX], fr.x[me[3][h1]][YY], fr.x[me[3][h1]][ZZ]);
                                                       int check9 =  angle6 (fr.x[me[2][h1]][XX], fr.x[me[2][h1]][YY], fr.x[me[2][h1]][ZZ], fr.x[me[3][h1]][XX], fr.x[me[3][h1]][YY], fr.x[mem[3][h1]][ZZ],
                                                                            fr.x[me[4][h1]][XX], fr.x[me[4][h1]][YY], fr.x[me[4][h1]][ZZ]);
                                                       int check10 =  angle6 (fr.x[me[3][h1]][XX], fr.x[me[3][h1]][YY], fr.x[me[3][h1]][ZZ], fr.x[me[4][h1]][XX], fr.x[me[4][h1]][YY], fr.x[mem[4][h1]][ZZ],
                                                                        fr.x[me[5][h1]][XX], fr.x[me[5][h1]][YY], fr.x[me[5][h1]][ZZ]);
                                                       int check11 =  angle6 (fr.x[me[4][h1]][XX], fr.x[me[4][h1]][YY], fr.x[me[4][h1]][ZZ], fr.x[me[5][h1]][XX], fr.x[me[5][h1]][YY], fr.x[mem[5][h1]][ZZ],
                                                                                                                               fr.x[me[0][h1]][XX], fr.x[me[0][h1]][YY], fr.x[me[0][h1]][ZZ]);
                                                        int check12 =  angle6 (fr.x[me[5][h1]][XX], fr.x[me[5][h1]][YY], fr.x[me[5][h1]][ZZ], fr.x[me[0][h1]][XX], fr.x[me[0][h1]][YY], fr.x[mem[0][h1]][ZZ],
                                                                                                                               fr.x[me[1][h1]][XX], fr.x[me[1][h1]][YY], fr.x[me[1][h1]][ZZ]);
                                                //   if (check7 == 1 && check8 == 1 && check9 == 1 && check10 == 1 && check11 == 1 && check12 == 1)
                                                                    {
                                                       float d2=sqrt(pow((fr.x[me[1][h1]][XX] - fr.x[me[5][h1]][XX]),2) + pow((fr.x[me[1][h1]][YY]- fr.x[me[5][h1]][YY]),2) + pow((fr.x[me[1][h1]][ZZ] - fr.x[me[5][h1]][ZZ]),2));
                                                       float d3=sqrt(pow((fr.x[me[2][h1]][XX] - fr.x[me[4][h1]][XX]),2) + pow((fr.x[me[2][h1]][YY]- fr.x[me[4][h1]][YY]),2) + pow((fr.x[me[2][h1]][ZZ] - fr.x[me[4][h1]][ZZ]),2));
                                                       float d4=sqrt(pow((fr.x[me[0][h1]][XX] - fr.x[me[3][h1]][XX]),2) + pow((fr.x[me[0][h1]][YY]- fr.x[me[3][h1]][YY]),2) + pow((fr.x[me[0][h1]][ZZ] - fr.x[me[3][h1]][ZZ]),2));
                                                         //    if (d2>0.35  || d2>0.35 || d3>0.35)
                                                                          {
                                                                            t3++;
                                                                          //  fprintf(fc," %d %d %d %d %d %d %d \n",me[0][h],me[1][h],me[2][h],me[3][h],me[4][h], me[5][h1],t3);           
                                                      A2[0][sixmemring] = me[0][h1];
                                                       A2[1][sixmemring] = me[1][h1];
                                                       A2[2][sixmemring] = me[2][h1];
                                                       A2[3][sixmemring] = me[3][h1];
                                                       A2[4][sixmemring] = me[4][h1];
                                                       A2[5][sixmemring] = me[5][h1];
                                                       A5[0][t3]=me[0][h1];
                                                       A5[1][t3]=me[1][h1];
                                                       A5[2][t3]=me[2][h1];
                                                       A5[3][t3]=me[3][h1];
                                                       A5[4][t3]=me[4][h1];
                                                       A5[5][t3]=me[5][h1];  
                                                       sixmemring++;
                                                    /*    index++;
                                                        index1++;

                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[me[0][h1]][XX]*10,fr.x[me[0][h1]][YY]*10,fr.x[me[0][h1]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[me[1][h1]][XX]*10,fr.x[me[1][h1]][YY]*10,fr.x[me[1][h1]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[me[2][h1]][XX]*10,fr.x[me[2][h1]][YY]*10,fr.x[me[2][h1]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[me[3][h1]][XX]*10,fr.x[me[3][h1]][YY]*10,fr.x[me[3][h1]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[me[4][h1]][XX]*10,fr.x[me[4][h1]][YY]*10,fr.x[me[4][h1]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[me[5][h1]][XX]*10,fr.x[me[5][h1]][YY]*10,fr.x[me[5][h1]][ZZ]*10);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-5,index);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-4,index-5);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-3,index-4);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-2,index-3);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-1,index-2);
                                                        fprintf(fq,"CONECT%5d%5d\n",index,index-1); */ 
                                                                                      
                                                                          }
                                                             }

                                                                                  }

                                                                         }
                                                                if (t3==2 && largecages==12   )
                                                                        {
                                                                           t7++;       
                                                                           f1++;
                                                                          gu8[f1]=gu2[largecages];                  
                                                           for (w6=1;w6<=largecages;w6++)
                                                                             {
                                                                             
                                                                  //     fprintf(fd," %d %d %d %d %d %d %d \n",A3[0][w6],A3[1][w6],A3[2][w6],A3[3][w6],A3[4][w6],w6,coun);
                                                                            index++;
                                                                            index1++;
                                                  /*    fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A3[0][w6]][XX]*10,fr.x[A3[0][w6]][YY]*10,fr.x[A3[0][w6]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A3[1][w6]][XX]*10,fr.x[A3[1][w6]][YY]*10,fr.x[A3[1][w6]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A3[2][w6]][XX]*10,fr.x[A3[2][w6]][YY]*10,fr.x[A3[2][w6]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A3[3][w6]][XX]*10,fr.x[A3[3][w6]][YY]*10,fr.x[A3[3][w6]][ZZ]*10);index++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A3[4][w6]][XX]*10,fr.x[A3[4][w6]][YY]*10,fr.x[A3[4][w6]][ZZ]*10);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-4,index);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-3,index-4);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-2,index-3);
                                                      fprintf(fq,"CONECT%5d%5d\n",index-1,index-2);
                                                      fprintf(fq,"CONECT%5d%5d\n",index,index-1); */
                                                                   for (w8=0;w8<=4;w8++)
                                                                               {
                                                                                n1++;
                                                                                  A4[n1]=A3[w8][w6];
                                                                                  gu3[n1]=gu2[w6];
                                                                  //  fprintf(fm," %d %d %d %d %d \n",A3[w8][w6],w6,coun,n1,gu2[w6]);    
                                                                               }   
                                                                            }                       
                                                           for (w7=1;w7<=t3;w7++)
                                                                             {                                                                     
                                                        index++;
                                                        index1++;               
                                                      /*  fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A5[0][w7]][XX]*10,fr.x[A5[0][w7]][YY]*10,fr.x[A5[0][w7]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A5[1][w7]][XX]*10,fr.x[A5[1][w7]][YY]*10,fr.x[A5[1][w7]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A5[2][w7]][XX]*10,fr.x[A5[2][w7]][YY]*10,fr.x[A5[2][w7]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A5[3][w7]][XX]*10,fr.x[A5[3][w7]][YY]*10,fr.x[A5[3][w7]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A5[4][w7]][XX]*10,fr.x[A5[4][w7]][YY]*10,fr.x[A5[4][w7]][ZZ]*10);index++;
                                                        fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index,index1,fr.x[A5[5][w7]][XX]*10,fr.x[A5[5][w7]][YY]*10,fr.x[A5[5][w7]][ZZ]*10);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-5,index);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-4,index-5);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-3,index-4);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-2,index-3);
                                                        fprintf(fq,"CONECT%5d%5d\n",index-1,index-2);
                                                        fprintf(fq,"CONECT%5d%5d\n",index,index-1);  */
                                                                        }
                                                                         }  
   /*     for (w9=1;w9<=n1;w9++)
                   {
                for (w10=w9+1;w10<=n1;)
                        {
                     if (A4[w9]==A4[w10])
                             {
                      for (w11=w10;w11<=n1;w11++)    
                                 {
                               A4[w11]=A4[w11+1];
                                 }
                               n1--;
                             } else                                
                               w10++;                       
                         }
                    } */
                      
                   for (w12=1;w12<=n1;w12++)   
                              {          
                                             
                    //     fprintf (fm, " %4d  \n",A4[w12]+1 );   
                    //      fprintf (fg, " %4d  \n",A4[w12]+2 ); 
                    //      fprintf (fg, " %4d  \n",A4[w12]+3 );
                              }
                        /*     if (gu3[n1]>0)
                                    {   
            //  fprintf (fa, " %4d \n",gu3[n1]+2 );                            
                                for (w13=0;w13<top.atoms.nr;w13++)
                                            {
                                         
              s1 = strcmp (*top.atoms.atomname[w13], "HW1");
              s2 = strcmp (*top.atoms.atomname[w13], "HW2");
                                                {
                                      if (s1==0 || s2==0)
                                                   { 
               distf=0;distg=0; 
               distf = sqrt(pow((fr.x[gu3[n1]+1][XX] - fr.x[w13][XX]),2) + pow((fr.x[gu3[n1]+1][YY] - fr.x[w13][YY]),2) + pow((fr.x[gu3[n1]+1][ZZ] - fr.x[w13][ZZ]),2));  
                distg = sqrt(pow((fr.x[gu3[n1]+2][XX] - fr.x[w13][XX]),2) + pow((fr.x[gu3[n1]+2][YY] - fr.x[w13][YY]),2) + pow((fr.x[gu3[n1]+2][ZZ] - fr.x[w13][ZZ]),2));
                                      if (distf<0.21 || distg<0.21)
                                                 {                                                 
                                                  u4++;  
                                                 trans[u4]=gu3[n1]+1;
                                                 trans1[u4]=gu3[n1]+2;
                                                 }
                                                  } 
                                                }
                                            }
                                    }   */                               
                           }
/*int qq=0,qq1=0;
      for (qq=1;qq<=u4;qq++)
              {
             int dd=0;
      for (qq1=1;qq1<=u4;qq1++)
              {
              if (qq!=qq1)
                  {
               if (trans[qq]==trans[qq1] || trans[qq]==trans1[qq1]|| trans1[qq]==trans1[qq1])      
                     {
                     dd++;
//  fprintf(fp," %f  %d %d %d \n",fr.time,trans[qq],trans1[qq],dd);                  
                     }
                  }
             }
              }*/
                                                                                                                                                       
int t5=0,in1[100000]={0};





         for(count10=0;count10<top.atoms.nr;count10++)     
                      {
      //    if (tim[count10]==fr.time)
                       {
  //           if (fr.x[count10][ZZ] >= 0 && fr.x[count10][ZZ] <= 6.16 || fr.x[count10][ZZ] >= 13.53 && fr.x[count10][ZZ] <= 19.47 )
                   {
             s1 = strcmp (*top.atoms.atomname[count10], "Ca");
             s2 = strcmp (*top.atoms.atomname[count10], "Ci");
               if (s2==0 || s1==0)
                  {
              w1++;
           in=0; 
                                for (count11=0; count11<top.atoms.nr; count11++)
                                                {
                                        s2 = strcmp (*top.atoms.atomname[count11], "OW1");
                                        s1 = strcmp (*top.atoms.atomname[count11], "OW");
                                           if (s2==0 || s1==0)
                                                  {
    
//  distMET = sqrt(pow(((x[count10]/10) - fr.x[count11][XX]),2) + pow(((y[count10]/10) - fr.x[count11][YY]),2) + pow(((z[count10]/10) - fr.x[count11][ZZ]),2));
//  distMET = sqrt(pow((fr.x[num[count10]][XX] - fr.x[count11][XX]),2) + pow((fr.x[num[count10]][YY] - fr.x[count11][YY]),2) + pow((fr.x[num[count10]][ZZ] - fr.x[count11][ZZ]),2));
    distMET = sqrt(pow((fr.x[count10][XX] - fr.x[count11][XX]),2) + pow((fr.x[count10][YY] - fr.x[count11][YY]),2) + pow((fr.x[count10][ZZ] - fr.x[count11][ZZ]),2));  
                                      if (distMET<=0.5890   )
                                                {
                                                        in++;
                                               ta2[w1][in]=count11;
                                               gu4[w1][in]=count10;
if (in==20)
{
inf1++;
for (coun=1; coun<=in; coun++)
                      {
                      ta3[inf1][coun]=ta2[w1][coun];           
                      gu5[inf1]=gu4[w1][coun]; 
                      in1[inf1]=in;  
                      }
}
                                                 }
                                                  }
                                                }
                  }
                     }
                  }
                   }
        for (coun=1; coun<=inf1; coun++)
                    { 
                      t1=0;t2=0;smallcages=0;n2=0; 
         for (count2=1; count2<=in1[coun]; count2++)
                      {
                             for (count3=1; count3<=in1[coun]; count3++)
                                                {
     distMET = sqrt(pow((fr.x[ta3[coun][count2]][XX]- fr.x[ta3[coun][count3]][XX]),2) + pow((fr.x[ta3[coun][count2]][YY] -fr.x[ta3[coun][count3]][YY]),2) + pow((fr.x[ta3[coun][count2]][ZZ] - fr.x[ta3[coun][count3]][ZZ]),2));
                                        if (distMET<=0.35 && distMET>0.25)
                                                {
                                                for (count5=1; count5<=in1[coun]; count5++)
                                                            {
         dista = sqrt(pow((fr.x[ta3[coun][count3]][XX] - fr.x[ta3[coun][count5]][XX]),2) + pow((fr.x[ta3[coun][count3]][YY]- fr.x[ta3[coun][count5]][YY]),2) + pow((fr.x[ta3[coun][count3]][ZZ] - fr.x[ta3[coun][count5]][ZZ]),2));
                                    if (dista<= 0.35 && dista>0.25) 
                                    
                                             if (ta3[coun][count2]!=ta3[coun][count5])
                                            {
                                                        for (count6=1; count6<=in1[coun]; count6++)
                                                            {
       distb = sqrt(pow((fr.x[ta3[coun][count5]][XX] - fr.x[ta3[coun][count6]][XX]),2) + pow((fr.x[ta3[coun][count5]][YY]- fr.x[ta3[coun][count6]][YY]),2) + pow((fr.x[ta3[coun][count5]][ZZ] - fr.x[ta3[coun][count6]][ZZ]),2));
                             if (distb<= 0.35 && distb>0.25 )
                             {
                                     if (ta3[coun][count3]!=ta3[coun][count6])
                                     {
                                                        for (count7=1; count7<=in1[coun]; count7++)
                                                            {
      distc = sqrt(pow((fr.x[ta3[coun][count6]][XX] - fr.x[ta3[coun][count7]][XX]),2) + pow((fr.x[ta3[coun][count6]][YY]- fr.x[ta3[coun][count7]][YY]),2) + pow((fr.x[ta3[coun][count6]][ZZ] - fr.x[ta3[coun][count7]][ZZ]),2));
                                  if (distc<= 0.35 && distc>0.25)
                                                   {
                                                          if (ta3[coun][count5]!=ta3[coun][count7])
                                                          {

                                                   for (count8=1; count8<=in1[coun]; count8++)
                                                            {
     distd = sqrt(pow((fr.x[ta3[coun][count7]][XX] - fr.x[ta3[coun][count8]][XX]),2) + pow((fr.x[ta3[coun][count7]][YY]- fr.x[ta3[coun][count8]][YY]),2) + pow((fr.x[ta3[coun][count7]][ZZ] - fr.x[ta3[coun][count8]][ZZ]),2));
                                     if (distd<= 0.35 && distd>0.25 )
                                                   {
                                                           if (ta3[coun][count6]!=ta3[coun][count8])
                                                           {
                                                        if (ta3[coun][count8]==ta3[coun][count2])
                                                        {
                                                           t1++; 
                                                                   inf++;
mem1[0][t1] = ta3[coun][count2]; mem1[1][t1] =ta3[coun][count3] ; mem1[2][t1] = ta3[coun][count5]; mem1[3][t1] =ta3[coun][count6]; mem1[4][t1] =ta3[coun][count7];
                                                        }
                                                                  for (count9=1; count9<=in1[coun]; count9++)
                                                                {
                       diste = sqrt(pow((fr.x[ta3[coun][count8]][XX] - fr.x[ta3[coun][count9]][XX]),2) + pow((fr.x[ta3[coun][count8]][YY]- fr.x[ta3[coun][count9]][YY]),2) + pow((fr.x[ta3[coun][count8]][ZZ] - fr.x[ta3[coun][count9]][ZZ]),2));

                                                                 if (diste<= 0.35 && diste>0.25  )
                                                                 {
                                                                      if (ta3[coun][count9]!=ta3[coun][count7])
                                                                      {

                                                                         if (ta3[coun][count9]==ta3[coun][count2])
                                                                         {
                                                                                 t2++;
                                                                                 inf2++;

                                        me1[0][t2] = ta3[coun][count2]; me1[1][t2] =ta3[coun][count3] ; me1[2][t2] = ta3[coun][count5]; me1[3][t2] =ta3[coun][count6]; me1[4][t2] =ta3[coun][count7]; me1[5][t2] =ta3[coun][count8];

                                                                         }
                                                                     }
                                                                  }
                                                                }
                                                        
                                                              }
                                                        }
                                                           }

                                                          }
                                                   }
                                                           }
                                      }
                             }
                                                           }
 
                                            }
                                    
                                                          }
                                                 }
                                             
                                                }
                                                }
                                                                    for(h=1;h<=t1;h++)
                                                                           {
                                                                                 check=0;
                                                                         for( t4=h+1;t4<=t1;t4++)
                                                                                  {
                                                                     counts=0;
                                                        for (e1 = 0; e1 <= 4; e1++)
                                                                {
                                                                         for (g = 0; g <=4; g++)
                                                                                  {
                                                                                if (mem1[e1][t4]==mem1[g][h])
                                                                                                {
                                                                                                    counts++;
                                                                                                }
                                                                                  }
                                                               if (counts >= 5)
                                                                      {
                                                                       check=-1;
                                                                      }
                                                                 }
                                                                       
                                                                                   }
                                                     int check_for6 = 0;
                                                         for (aj = 1; aj <= t2; aj++)
                                                          {
                                                           count11 = 0;
                                                           for (ac = 0; ac <= 5; ac++)
                                                            {
                                                             for (vc = 0; vc <= 4; vc++)
                                                             {
                                                              if (me1[ac][aj] == mem1[vc][h])
                                                               {
                                                                 count11++;
                                                               }
                                                             }
                                                             if (count11 >= 5)
                                                             {
                                                              check_for6 = -1;
                                                             }
                                                            }
                                                           }
                                                                   if(check == 0 && check_for6==0  )
                                                                       {
                                                       int check1 =  angle5 (fr.x[mem1[0][h]][XX], fr.x[mem1[0][h]][YY], fr.x[mem1[0][h]][ZZ], fr.x[mem1[1][h]][XX], fr.x[mem1[1][h]][YY], fr.x[mem1[1][h]][ZZ],
                                                              fr.x[mem1[2][h]][XX], fr.x[mem1[2][h]][YY], fr.x[mem1[2][h]][ZZ]);
                                                       int check2 =  angle5 (fr.x[mem1[1][h]][XX], fr.x[mem1[1][h]][YY], fr.x[mem1[1][h]][ZZ], fr.x[mem1[2][h]][XX], fr.x[mem1[2][h]][YY], fr.x[mem1[2][h]][ZZ],
                                                                                                                     fr.x[mem1[3][h]][XX], fr.x[mem1[3][h]][YY], fr.x[mem1[3][h]][ZZ]);
                                                       int check3 =  angle5 (fr.x[mem1[2][h]][XX], fr.x[mem1[2][h]][YY], fr.x[mem1[2][h]][ZZ], fr.x[mem1[3][h]][XX], fr.x[mem1[3][h]][YY], fr.x[mem1[3][h]][ZZ],
                                                                            fr.x[mem1[4][h]][XX], fr.x[mem1[4][h]][YY], fr.x[mem1[4][h]][ZZ]);
                                                       int check4 =  angle5 (fr.x[mem1[3][h]][XX], fr.x[mem1[3][h]][YY], fr.x[mem1[3][h]][ZZ], fr.x[mem1[4][h]][XX], fr.x[mem1[4][h]][YY], fr.x[mem1[4][h]][ZZ],
                                                                        fr.x[mem1[0][h]][XX], fr.x[mem1[0][h]][YY], fr.x[mem1[0][h]][ZZ]);
                                                       int check5 =  angle5 (fr.x[mem1[4][h]][XX], fr.x[mem1[4][h]][YY], fr.x[mem1[4][h]][ZZ], fr.x[mem1[0][h]][XX], fr.x[mem1[0][h]][YY], fr.x[mem1[0][h]][ZZ],
                                                                                                                               fr.x[mem1[1][h]][XX], fr.x[mem1[1][h]][YY], fr.x[mem1[1][h]][ZZ]);
                                                //  if (check1 == 1 && check2 == 1 && check3 == 1 && check4 == 1 && check5 == 1 )

                                                                  {
                                                       float d1=sqrt(pow((fr.x[mem1[1][h]][XX] - fr.x[mem1[4][h]][XX]),2) + pow((fr.x[mem1[1][h]][YY]- fr.x[mem1[4][h]][YY]),2) + pow((fr.x[mem1[1][h]][ZZ] - fr.x[mem1[4][h]][ZZ]),2));
                                                       float d4=sqrt(pow((fr.x[mem1[2][h]][XX] - fr.x[mem1[4][h]][XX]),2) + pow((fr.x[mem1[2][h]][YY]- fr.x[mem1[4][h]][YY]),2) + pow((fr.x[mem1[2][h]][ZZ] - fr.x[mem1[4][h]][ZZ]),2));
                                                       float d5=sqrt(pow((fr.x[mem1[1][h]][XX] - fr.x[mem1[3][h]][XX]),2) + pow((fr.x[mem1[1][h]][YY]- fr.x[mem1[3][h]][YY]),2) + pow((fr.x[mem1[1][h]][ZZ] - fr.x[mem1[3][h]][ZZ]),2));
                                                       float d6=sqrt(pow((fr.x[mem1[0][h]][XX] - fr.x[mem1[3][h]][XX]),2) + pow((fr.x[mem1[0][h]][YY]- fr.x[mem1[3][h]][YY]),2) + pow((fr.x[mem1[0][h]][ZZ] - fr.x[mem1[3][h]][ZZ]),2));  
                                                       float d7=sqrt(pow((fr.x[mem1[0][h]][XX] - fr.x[mem1[2][h]][XX]),2) + pow((fr.x[mem1[0][h]][YY]- fr.x[mem1[2][h]][YY]),2) + pow((fr.x[mem1[0][h]][ZZ] - fr.x[mem1[2][h]][ZZ]),2));
                                                             //      if (d1>0.35 || d4>0.35 || d5>0.35 || d6>0.35 || d7>0.35)      
                                                                           {          
                                                                       smallcages++;
                                                                        //  fprintf(fc," %d %d %d %d %d %d %d \n",mem[0][h],mem[1][h],mem[2][h],mem[3][h],mem[4][h],smallcages,coun);
                                                       A1[0][fivememring] = mem1[0][h];
                                                       A1[1][fivememring] = mem1[1][h];
                                                       A1[2][fivememring] = mem1[2][h];
                                                       A1[3][fivememring] = mem1[3][h];
                                                       A1[4][fivememring] = mem1[4][h];
                                                       A6[0][smallcages]=mem1[0][h];
                                                       A6[1][smallcages]=mem1[1][h];
                                                       A6[2][smallcages]=mem1[2][h];
                                                       A6[3][smallcages]=mem1[3][h];
                                                       A6[4][smallcages]=mem1[4][h];
                                                       gu6[smallcages]=gu5[coun];
                                                       fivememring++;
                                                       index2++;
                                                       index3++;
                                           /*           fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[mem[0][h]][XX]*10,fr.x[mem[0][h]][YY]*10,fr.x[mem[0][h]][ZZ]*10);index2++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[mem[1][h]][XX]*10,fr.x[mem[1][h]][YY]*10,fr.x[mem[1][h]][ZZ]*10);index2++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[mem[2][h]][XX]*10,fr.x[mem[2][h]][YY]*10,fr.x[mem[2][h]][ZZ]*10);index2++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[mem[3][h]][XX]*10,fr.x[mem[3][h]][YY]*10,fr.x[mem[3][h]][ZZ]*10);index2++;
                                                      fprintf(fq,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[mem[4][h]][XX]*10,fr.x[mem[4][h]][YY]*10,fr.x[mem[4][h]][ZZ]*10);
                                                      fprintf(fq,"CONECT%5d%5d\n",index2-4,index2);
                                                      fprintf(fq,"CONECT%5d%5d\n",index2-3,index2-4);
                                                      fprintf(fq,"CONECT%5d%5d\n",index2-2,index2-3);
                                                      fprintf(fq,"CONECT%5d%5d\n",index2-1,index2-2);
                                                      fprintf(fq,"CONECT%5d%5d\n",index2,index2-1);    */                                                                               
                                                                           }
                                                                  }
                                                              }                        
                                                                 }
                                                        if ( smallcages==12     )
                                                            {
                                                            t5++;
                                                        for (w14=1;w14<=smallcages;w14++)
                                                                {
                                                                            index2++;
                                                                            index3++;
                                                 /*     fprintf(fb,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[A6[0][w14]][XX]*10,fr.x[A6[0][w14]][YY]*10,fr.x[A6[0][w14]][ZZ]*10);index2++;
                                                      fprintf(fb,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[A6[1][w14]][XX]*10,fr.x[A6[1][w14]][YY]*10,fr.x[A6[1][w14]][ZZ]*10);index2++;
                                                      fprintf(fb,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[A6[2][w14]][XX]*10,fr.x[A6[2][w14]][YY]*10,fr.x[A6[2][w14]][ZZ]*10);index2++;
                                                      fprintf(fb,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[A6[3][w14]][XX]*10,fr.x[A6[3][w14]][YY]*10,fr.x[A6[3][w14]][ZZ]*10);index2++;
                                                      fprintf(fb,"ATOM  %5d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",index2,index3,fr.x[A6[4][w14]][XX]*10,fr.x[A6[4][w14]][YY]*10,fr.x[A6[4][w14]][ZZ]*10);
                                                      fprintf(fb,"CONECT%5d%5d\n",index2-4,index2);
                                                      fprintf(fb,"CONECT%5d%5d\n",index2-3,index2-4);
                                                      fprintf(fb,"CONECT%5d%5d\n",index2-2,index2-3);
                                                      fprintf(fb,"CONECT%5d%5d\n",index2-1,index2-2);
                                                      fprintf(fb,"CONECT%5d%5d\n",index2,index2-1);  */
                                                                   for (w15=0;w15<=4;w15++)
                                                                               {
                                                                                n2++;
                                                                                  A7[n2]=A6[w15][w14];
                                                                                  gu7[n2]=gu6[w14];
                                                                               }                                                                                                                                 
                                                                }
                                                                           
                                                            }
      /*  for (w16=1;w16<=n2;w16++)
                   {
                for (w17=w16+1;w17<=n2;)
                        {
                     if (A7[w16]==A7[w17])
                             {
                      for (w18=w17;w18<=n2;w18++)
                                 {
                               A7[w18]=A7[w18+1];
                                 }
                               n2--;
                             } else
                               w17++;
                         }
                    }*/
                   for (w19=1;w19<=n2;w19++)
                              {

                     //    fprintf (fm, " %4d  \n",A7[w19]+1 );    
                     //     fprintf (fg, " %4d  \n",A7[w19]+2 );  
                     //     fprintf (fg, " %4d  \n",A7[w19]+3 );
                              }
                          /*   if (gu7[n2]>0)
                                    {
           //   fprintf (fa, " %4d \n",gu7[n2]+2 );                             
                                for (w20=0;w20<top.atoms.nr;w20++)
                                            {

              s1 = strcmp (*top.atoms.atomname[w20], "HW1");
              s2 = strcmp (*top.atoms.atomname[w20], "HW2");
                                                {
                                      if (s1==0 || s2==0)
                                                   {
               distf=0;distg=0;
               distf = sqrt(pow((fr.x[gu7[n2]+1][XX] - fr.x[w20][XX]),2) + pow((fr.x[gu7[n2]+1][YY] - fr.x[w20][YY]),2) + pow((fr.x[gu7[n2]+1][ZZ] - fr.x[w20][ZZ]),2));
                distg = sqrt(pow((fr.x[gu7[n2]+2][XX] - fr.x[w20][XX]),2) + pow((fr.x[gu7[n2]+2][YY] - fr.x[w20][YY]),2) + pow((fr.x[gu7[n2]+2][ZZ] - fr.x[w20][ZZ]),2));
                                      if (distf<0.21 || distg<0.21)
                                                 {
                                                  u5++;
                                                 }
                                                  }
                                                }
                                            }
                                    }*/                                                                          
                           }        

//fprintf (fq, "TER \n");
//fprintf (fq, "ENDMDL \n");
                                                                                                                                                                                                                   
fprintf(fd," %f %d \n",fr.time,t7+t5);
//fprintf(fp," %f  %d %d \n",fr.time,u5,u4);

}while(read_next_frame(status,&fr));



  thanx(stderr);
  
  return 0;
}

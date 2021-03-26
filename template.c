static char *SRCID_template_c = "$Id$";
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "this is a small test program meant to serve as a template ",
    "select which atom you want to examine with the -n argument."
  };
  static int n=1;
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    }
  };
  t_topology top;
  int        ePBC;
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
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  sfree(xtop);
  n=n-1; /* Our enumeration started on 1, but C starts from 0 */
  /* check that this atom exists */
  if(n<0 || n>(top.atoms.nr)) 
  {
    printf("Error: Atom number %d is out of range.\n",n);
    exit(1);
  }
  printf("Atom name: %s\n",*(top.atoms.atomname[n]));
  printf("Atom charge: %f\n",top.atoms.atom[n].q);
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
  FILE *out;
  out=fopen("fromed.pdb","wb");
  FILE *out1;
  out1=fopen("mcg.txt","wb");
  FILE *out2;
  out2=fopen("mcgOP.txt","wb");
  FILE *out3;
  out3=fopen("1large.txt","wb");
  FILE *out4;
  out4=fopen("2large.txt","wb");
  FILE *out5;
  out5=fopen("ones.txt","wb");
  FILE *out6;
  out6=fopen("two.txt","wb");
  FILE *out7;
  out7=fopen("three.txt","wb");
  FILE *out8;
  out8=fopen("four.txt","wb");
  FILE *out9;
  out9=fopen("five.txt","wb");

int zai=1;
/*int count=0;  */
do {
fprintf(out,"%f\n",fr.time);
int  group[4][5000],NG=0, group2[4][5000],weight[2][1000],comb[2][10000],combt[2][10000],cb=0,wt=0;
int i=0,k=0,s1=-1,s2=-1,s3=-1,s4=-1,s5=-1,s6=-1;float cosphi=0;int n=0,count=0,p1=0,p2=0,p3=0,p4=0,p5=0,p6=0,p7=0,p8=0,p9=0,p10=0,knum=0;
while (i<top.atoms.nr)
{
 s1=strcmp (*top.atoms.atomname[i],"CT");
 s2=strcmp (*top.atoms.atomname[i],"Cl");
if (s1==0 ||s2==0)
{
int k=i+1;
while (k<top.atoms.nr)
{
if (k != i)
{
 s3=strcmp (*top.atoms.atomname[k],"Cl");s4=strcmp (*top.atoms.atomname[k],"CT");
float d=0;
if (s3==0 || s4==0)
{
float x=fr.x[i][XX]-fr.x[k][XX];float y=fr.x[i][YY]-fr.x[k][YY];float z=fr.x[i][ZZ]-fr.x[k][ZZ];
d=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
int Wnum=0;
if (d<=0.87 && d>=0.5)
{
int j=0,Wa[1][6];
int Waa=0,Wbb=0,Wcc=0,Wdd=0,Wee=0;float Paa=0,Pbb=0,Pcc=0,Pdd=0,Pee=0,Qaa=0,Qbb=0,Qcc=0,Qdd=0,Qee=0;
while(j<top.atoms.nr) 
{
 s5=strcmp (*top.atoms.atomname[j],"OW3"); s6=strcmp (*top.atoms.atomname[j],"OW");
float xmid=0,ymid=0,zmid=0,Rcut=0.45,Rco=0,dxr=0,dyr=0,dzr=0,Angmin=0.50,Angmax=1;
if (s5==0 || s6==0)
{
float xmid=(fr.x[i][XX]+fr.x[k][XX])/2;float ymid=(fr.x[i][YY]+fr.x[k][YY])/2;float zmid=(fr.x[i][ZZ]+fr.x[k][ZZ])/2;
float dxr=fr.x[j][XX]-xmid;float dyr=fr.x[j][YY]-ymid;float dzr=fr.x[j][ZZ]-zmid;
Rco=sqrt(pow(dxr,2)+pow(dyr,2)+pow(dzr,2));
 
float Tx=0,Ty=0,Tz=0,Px=0,Py=0,Pz=0,Qx=0,Qy=0,Qz=0,Rx=0,Ry=0,Rz=0,TT=0,PP=0,QQ=0,RR=0,TP=0,RQ=0,TTPP=0,RRQQ=0,ang1=0,ang2=0,TcP=0,RcQ=0,Amax=0.93969,fmin=0.5,gmin=0.70710678118,Amin=0.57357643;
if (Rco<=Rcut)
{
float Tx=(fr.x[k][XX]-fr.x[i][XX]);float Ty=(fr.x[k][YY]-fr.x[i][YY]);float Tz=(fr.x[k][ZZ]-fr.x[i][ZZ]);
float Px=(fr.x[j][XX]-fr.x[i][XX]);float Py=(fr.x[j][YY]-fr.x[i][YY]);float Pz=(fr.x[j][ZZ]-fr.x[i][ZZ]);
float Qx=(fr.x[j][XX]-fr.x[k][XX]);float Qy=(fr.x[j][YY]-fr.x[k][YY]);float Qz=(fr.x[j][ZZ]-fr.x[k][ZZ]);float Rx=-1*Tx;float Ry=-1*Ty;float Rz=-1*Tz;
float TT=sqrt(pow(Tx,2)+pow(Ty,2)+pow(Tz,2));float PP=sqrt(pow(Px,2)+pow(Py,2)+pow(Pz,2));
float RR=sqrt(pow(Rx,2)+pow(Ry,2)+pow(Rz,2));float QQ=sqrt(pow(Qx,2)+pow(Qy,2)+pow(Qz,2));
float TP=Tx*Px+Ty*Py+Tz*Pz;float RQ=Rx*Qx+Ry*Qy+Rz*Qz;float ang1=TP/(TT*PP);float ang2=RQ/(RR*QQ);
if(ang1<=Amax && ang1>=Amin)
{
if(ang2<=Amax && ang2>=Amin)
{
if (Wnum==0){ Paa=ang1;Qaa=ang2;} if (Wnum==1){ Pbb=ang1;Qbb=ang2;} if (Wnum==2){ Pcc=ang1;Qcc=ang2;}if (Wnum==3){ Pdd=ang1;Qdd=ang2;}if (Wnum==4){Pee=ang1;Qee=ang2;}
/*if (Wnum==0){ Waa=j;} if (Wnum==1){ Wbb=j;} if (Wnum==2){ Wcc=j;}if (Wnum==3){ Wdd=j;}if (Wnum==4){ Wee=j;}*/
Wnum++;
if (Wnum>=5) 
{
if(fr.time==1500.000)
{
fprintf(out8," %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n",Paa,Qaa,Pbb,Qbb,Pcc,Qcc,Pdd,Qdd,Pee,Qee);
}
/* second searchAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA */
/*fprintf(out, " %f  %d %d %d \n",fr.time,i,k,cb);*/
comb[0][cb]=i;comb[1][cb]=k;
combt[0][cb]=i;combt[1][cb]=k;
/*if (fr.time==0)
{
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",i+1,*top.atoms.atomname[i],*top.atoms.resname[top.atoms.atom[i].resnr],top.atoms.atom[i].resnr+1,(fr.x[i][XX])*10,(fr.x[i][YY])*10,(fr.x[i][ZZ])*10);
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",k+1,*top.atoms.atomname[k],*top.atoms.resname[top.atoms.atom[k].resnr],top.atoms.atom[k].resnr+1,(fr.x[k][XX])*10,(fr.x[k][YY])*10,(fr.x[k][ZZ])*10);
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",Waa+1,*top.atoms.atomname[Waa],*top.atoms.resname[top.atoms.atom[Waa].resnr],top.atoms.atom[Waa].resnr+1,(fr.x[Waa][XX])*10,(fr.x[Waa][YY])*10,(fr.x[Waa][ZZ])*10);
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",Wbb+1,*top.atoms.atomname[Wbb],*top.atoms.resname[top.atoms.atom[Wbb].resnr],top.atoms.atom[Wbb].resnr+1,(fr.x[Wbb][XX])*10,(fr.x[Wbb][YY])*10,(fr.x[Wbb][ZZ])*10);
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",Wcc+1,*top.atoms.atomname[Wcc],*top.atoms.resname[top.atoms.atom[Wcc].resnr],top.atoms.atom[Wcc].resnr+1,(fr.x[Wcc][XX])*10,(fr.x[Wcc][YY])*10,(fr.x[Wcc][ZZ])*10);
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",Wdd+1,*top.atoms.atomname[Wdd],*top.atoms.resname[top.atoms.atom[Wdd].resnr],top.atoms.atom[Wdd].resnr+1,(fr.x[Wdd][XX])*10,(fr.x[Wdd][YY])*10,(fr.x[Wdd][ZZ])*10);
fprintf(out2,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",Waa+1,*top.atoms.atomname[Wee],*top.atoms.resname[top.atoms.atom[Waa].resnr],top.atoms.atom[Wee].resnr+1,(fr.x[Wee][XX])*10,(fr.x[Wee][YY])*10,(fr.x[Wee][ZZ])*10);
}
*/

cb++;
knum++;
/* second search AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*/
}}}}}j++;}}}}k++;}
if (knum>=1)
{

/*fprintf(out2, " %f  %d %d %d  \n",fr.time,i,knum,wt);*/
/*weight[0][wt]=i;weight[1][wt]=knum;*/
wt++;knum=0;
}
}i++;
}
count=count-1;

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
int mapsize=2*cb+10;
int lsize=mapsize/2;
int map[mapsize][mapsize];
int Csize[2][mapsize];
int jadoo=0;
int appy=0;
for(int mk=0;mk<=mapsize;mk++)
{
for(int asj=0;asj<=lsize;asj++)
{
if(comb[0][asj]!=0||comb[1][asj]!=0) 
{
map[0][mk]=comb[0][asj];
map[1][mk]=comb[1][asj];
comb[0][asj]=0;
comb[1][asj]=0;
asj=lsize;
}
 }
int si=2,sj=3;
for(int mj=0;mj<=mapsize;mj++)
{
for(int nj=0;nj<=lsize;nj++)
{
if(comb[0][nj]!=0||comb[1][nj]!=0)
{
if(comb[0][nj]==map[mj][mk]||comb[1][nj]==map[mj][mk])
{
map[si][mk]=comb[0][nj];
map[sj][mk]=comb[1][nj];
/*printf("%f %d %d %d %d %d\n",fr.time,jadoo,si,mk,comb[0][nj],comb[1][nj]);*/
comb[0][nj]=0;
comb[1][nj]=0;
si=si+2;
sj=sj+2;
}
}
}
if(map[mj][mk]==0)
{
mj=mapsize;
}
}
Csize[0][mk]=mk;
Csize[1][mk]=sj;
jadoo=mk;
/*printf("%f %d %d\n",fr.time,Csize[1][mk],Csize[0][mk]);*/
/*printf("%f %d %d %d\n",fr.time,jadoo,si,mk);*/
if(map[1][mk]==0)
{
mk=mapsize;
}
}
fprintf(out5,"%f  %d\n",fr.time,jadoo+1);
/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
/*finding maxmum arrysize*/
int gi=0,gj=0,secj=0,seci=0;
for(int xi=0;xi<=jadoo;xi++)
{
if( gi<Csize[1][xi])
{
seci=gi;
gi=Csize[1][xi];
secj=gj;
gj=Csize[0][xi];
}
}
gi=gi;

/*int pi=0,pj=0;
for(int xi=0;xi<=jadoo;xi++)
{
if(Csize[1][xi]<gi)
{
if( pi<Csize[1][xi])
{
pi=Csize[1][xi];
pj=Csize[0][xi];*/
/*printf("%f %d %d %d\n",fr.time,Csize[1][xi],Csize[0][xi],jadoo);*/
/*}
}
}
fprintf( out4," %f  %d  \n",fr.time,pi);*/
/*********************/
int ei=0,ej=0;
for(int xi=0;xi<=jadoo;xi++)
{


ei=Csize[1][xi];
ej=Csize[0][xi];
fprintf(out4, " %d    %d\n", ei, zai);

}

/********************/

/*printf("%f %d %d %d\n",fr.time,gi,jadoo,gj);*/
/*&&&&&&&&&&&&&&&&&&&&&*/
/*sorting of maximum array*/
int happy=0;
/*fprintf(out, "REMARK    GENERATED BY TRJCONV\nTITLE     co2_hydrate t=   %f\nCRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P %d\nMODEL %8d\n",fr.time,fr.box[XX][XX]*10,fr.box[YY][YY]*10,fr.box[ZZ][ZZ]*10,90.00,90.00,90.00,1,zai);*/

for(int xj=0;xj<=gi;xj++)
{
int big=0;
for(int xk=0;xk<=gi;xk++)
{
if(big<map[xk][gj])
{
big=map[xk][gj];
}
}
for(int xk=0;xk<=gi;xk++)
{
if(big==map[xk][gj])
{
map[xk][gj]=0;
}
}
if(big!=0)
{
/*fprintf(out,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",big+1,*top.atoms.atomname[big],*top.atoms.resname[top.atoms.atom[big].resnr],top.atoms.atom[big].resnr+1,(fr.x[big][XX])*10,(fr.x[big][YY])*10,(fr.x[big][ZZ])*10);*/
happy=happy+1;
}
if(big==0)
{
fprintf(out, "TER\nENDMDL\n");
xj=gi;
}
}
int happysec=0;
if (secj!=gj)
{
for(int xp=0;xp<=seci;xp++)
{
int bigsec=0;
for(int xq=0;xq<=seci;xq++)
{
if(bigsec<map[xq][secj])
{
bigsec=map[xq][secj];
}
}
for(int xq=0;xq<=seci;xq++)
{
if(bigsec==map[xq][secj])
{
map[xq][secj]=0;
}
}
if(bigsec!=0)
{
/*fprintf(out,"ATOM  %5d  %-3s %-3s%6d    %8.3f%8.3f%8.3f  1.00  0.00\n",big+1,*top.atoms.atomname[big],*top.atoms.resname[top.atoms.atom[big].resnr],top.atoms.atom[big].resnr+1,(fr.x[big][XX])*10,(fr.x[big][YY])*10,(fr.x[big][ZZ])*10);*/
happysec=happysec+1;
}
if(bigsec==0)
{
fprintf(out, "TER\nENDMDL\n");
xp=seci;
}
}
}

int humpty=0;
/*total methane in all cluster**/
for (int hk=0;hk<=lsize;hk++)
{
int bingo=0;
for (int hi=0;hi<=lsize;hi++)
{
for(int hj=0;hj<=1;hj++)
{
if(combt[hj][hi]!=0)
{
if (bingo<combt[hj][hi])
{
bingo=combt[hj][hi];
}}}}
humpty=humpty+1;
if (bingo==0)
{
hk=lsize;
}
for(int di=0;di<=lsize;di++)
{
for(int dj=0;dj<=1;dj++)
{

if(combt[dj][di]==bingo)
{
combt[dj][di]=0;
}}
}} 

/*********************/




fprintf(out7, " %f  %d  \n",fr.time,humpty);

zai=zai+1; 
fprintf(out6, " %f  %d  \n",fr.time,happy);
fprintf(out9, " %f  %d  \n",fr.time,happysec);
/*&&&&&&&&&&&&&&&&&&&&&&&*/
  } while(read_next_frame(status,&fr));
/*fclose(out1);*/
  thanx(stderr);
  
  return 0;
}


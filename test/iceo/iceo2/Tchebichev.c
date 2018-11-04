

double Chebychev(double *tmp, int D)                 /*  Chebychev Polynomial */
{                                                    /* Valid for D=9 or D=17 */
   int i,j;
   double px,x=-1,result=0,dx;
   
   if ( D == 9 )
     dx=72.66066;
   if ( D == 17 )
     dx=10558.145;

   for (i=0;i<=100;i++)
   {
      px=tmp[0];
      for (j=1;j<D;j++) px=x*px+tmp[j];
      if (px<-1 || px>1) result+=(1.-px)*(1.-px);
      x+=.02;
   }
   px=tmp[0];
   for (j=1;j<D;j++) px=1.2*px+tmp[j];
   px=px-dx;
   if (px<0) result+=px*px;
   px=tmp[0];
   for (j=1;j<D;j++) px=-1.2*px+tmp[j];
   px=px-dx;
   if (px<0) result+=px*px;
   return result;
}








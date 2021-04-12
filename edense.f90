module edense
implicit none
double precision,parameter::hbarcc=197.32697, pii=3.1415926535897932384626433, mee=0.511d0,en=2.7182818284590452


contains


double precision function edens(p,temp,xa)
        
        double precision::p,temp,z,x,triu,integrate,y,dif,dif1
        double precision,dimension(3000):: denomi1,denomi2,numfac,xx
        integer::i,j
        logical::xa

        x=1d-200
        do i=1,3000
        denomi1(i)=x**2
        denomi2(i)=(x**4+1.)/2.
        numfac(i)=-4*temp**2/hbarcc**3/pii**2*x*log(x)*sqrt(4*log(x)**2*temp**2-mee**2)
        x=x+exp(-0.5*mee/temp)/3000
        end do
        
        x=0d0
        dif=1d99
        if (xa) then

        do while(dif>=0)
        triu=0d0
                do j=1,3000
                triu=triu+(exp(x)-exp(-x))/2*numfac(j)/(denomi1(j)*(exp(x)+exp(-x))/2+denomi2(j))*exp(-mee*0.5/temp)/3001 
                end do    
        dif=p-triu
        x=x+0.01
        end do

                do while(dif<=0)
                triu=0d0
                        do j=1,3000
                        triu=triu+(exp(x)-exp(-x))/2*numfac(j)/(denomi1(j)*(exp(x)+exp(-x))/2+denomi2(j))*exp(-mee*0.5/temp)/3001 
                        end do    
                        dif=p-triu
                        x=x-0.001
                end do
        
        z=exp(x)
        else
        z=1d0
        end if        
        integrate=0d0
        y=1d-200
        do i=1,3000
        integrate=integrate+numfac(i)*1/(1/z+denomi1(i))*exp(-mee*0.5/temp)/3000
        y=y+exp(-mee*0.5/temp)/3000        
        end do

        edens=integrate     
end function edens

end module edense

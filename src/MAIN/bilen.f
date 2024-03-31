      subroutine bilen(ax,ay,aa,x,y,a,dadx,dady,lx,ly,i,j)
c..bilinear interpolation with derivatives
c..   ax = vector for x axis
c..   ay = vector for y axis
c..   aa = 2D matrix of function values
c..   x  = x variable value
c..   y  = y variable value
c..   a  = interpolated function value
c..   dadx = derivitive of interpolant with respect to x
c..   dady = derivitive of interpolant with respect to y
c..   lx   = x dimension of aa array
c..   ly   = y dimension of aa array
c..   i    = nearest lower index in x
c..   j     = nearest lower index in y
      implicit none
      integer*4 lx, ly
      integer*4 i, j
      real*8 aa(lx,ly), ax(lx), ay(ly)
      real*8 x, y, a, t, u
      real*8 dtdx, dudy, dadx, dady
c---------------------------------------------------------------------
c..interpolation factors
      dtdx = 1.0d0/( ax(i+1) - ax(i) )
      t    = (x - ax(i) )*dtdx

      dudy = 1.0d0/( ay(j+1) - ay(j) )
      u    = (y - ay(j) )*dudy

c..bilinear interpolation function
      a = (1.0d0 - t)*(1.0d0 - u)*aa(i  ,j)
     1  +          t *(1.0d0 - u)*aa(i+1,j)
     2  +          t *         u *aa(i+1,j+1)
     3  + (1.0d0 - t)*         u *aa(i  ,j+1)

c..analytic derivatives of approximation function a
      dadx = dtdx*( (1.0d0 - u)*( -aa(  i,j  ) + aa(i+1,j  ) )
     1             +         u *(  aa(i+1,j+1) - aa(i  ,j+1) ) )
      dady = dudy*( (1.0d0 - t)*( -aa(  i,j  ) + aa(i  ,j+1) )
     1             +         t *( -aa(i+1,j  ) + aa(i+1,j+1) ) )

      return

      end


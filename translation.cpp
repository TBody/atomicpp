struct return_regrid{
      int nx;
      int ny;
      vector<double> tx;
      vector<double> ty;
      vector<double> c;
      double fp;
      int ier;
};

return_regrid regrid_smth(std::vector<double> x, std::vector<double> y, std::vector<std::vector<double>> z){
      
      // Translation of interpolate/src/fitpack.pyf
      // nx,tx,ny,ty,c,fp,ier = regrid_smth(x,y,z,[xb,xe,yb,ye,kx,ky,s])

      int iopt = 0;
      int mx = x.size();
      int my = y.size();

      int kx = 3;
      int ky = 3;  

      if (not((mx > kx) and (my > ky))){
            throw runtime_error("Grid too small for bicubic interpolation")
      };

      double s = 0.0;    

      double xb = std::min(x,mx);
      double xe = std::max(x,mx);
      double yb = std::min(y,my);
      double ye = std::max(y,my);      

      int nxest = mx+kx+1;
      int nyest = my+ky+1;

      std::vector<double> tx(nxest, 0.0);
      std::vector<double> ty(nyest, 0.0);
      std::vector<double> c((nxest-kx-1)*(nyest-ky-1), 0.0);
      double fp = 0.0;

      int lwrk = 4 + nxest * (my + 2 * kx + 5) + nyest * (2 * ky + 5) + mx * (kx + 1) + my * (ky + 1) + max(my, nxest);
      std::vector<double> wrk(lwrk, 0.0);
      int kwrk = 3 + mx + my + nxest + nyest;
      std::vector<int> iwrk(kwrk, 0);

      int ier = 0;
      // End of translation of interpolate/src/fitpack.pyf

}

















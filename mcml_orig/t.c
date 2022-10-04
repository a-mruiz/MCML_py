void ScaleA(InputStruct In_Parm, OutStruct * Out_Ptr)
{
  short nz = In_Parm.nz;
  short nr = In_Parm.nr;
  double dz = In_Parm.dz;
  double dr = In_Parm.dr;
  short nl = In_Parm.num_layers;
  short iz,ir;
  short il;
  double scale1;
  
  /* Scale A_rz. */
  scale1 = 2.0*PI*dr*dr*dz*In_Parm.num_photons;	
	/* volume is 2*pi*(ir+0.5)*dr*dr*dz.*/ 
	/* ir+0.5 to be added. */
  for(iz=0; iz<nz; iz++) 
    for(ir=0; ir<nr; ir++) 
      Out_Ptr->A_rz[ir][iz] /= (ir+0.5)*scale1;
  
  /* Scale A_z. */
  scale1 = 1.0/(dz*In_Parm.num_photons);
  for(iz=0; iz<nz; iz++) 
    Out_Ptr->A_z[iz] *= scale1;
  
  /* Scale A_l. Avoid int/int. */
  scale1 = 1.0/(double)In_Parm.num_photons;	
  for(il=0; il<=nl+1; il++)
    Out_Ptr->A_l[il] *= scale1;
  
  Out_Ptr->A *=scale1;
}
/***********************************************************
 *	Return the index to the layer according to the index
 *	to the grid line system in z direction (Iz).
 *
 *	Use the center of box.
 ****/
short IzToLayer(short Iz, InputStruct In_Parm)
{
  short i=1;	/* index to layer. */
  short num_layers = In_Parm.num_layers;
  double dz = In_Parm.dz;
  
  while( (Iz+0.5)*dz >= In_Parm.layerspecs[i].z1 
	&& i<num_layers) i++;
  
  return(i);
}
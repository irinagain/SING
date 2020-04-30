makecifti <- function(table, filename, template = '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii', toolboxloc1 = '~/Dropbox/Applications2/hcp_ciftimatlab', toolboxloc2 = '~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8',wbloc = '~/Applications2/workbench/bin_macosx64/wb_command',scratchloc = '~',matlabloc = '/Applications/MATLAB_R2018b.app/bin/matlab'){
  
  # Adapted from https://mandymejia.wordpress.com/2016/10/28/r-function-to-write-cifti-files/ 
  # table = Vxp matrix
  # table can be a vector if only a single CIFTI is needed
  # V = number of locations in both or one hemisphere
  # p = number of time points in .dtseries.nii
  # filename = vector length p, containing file names to be written (no extension)
  # template = name of existing CIFTI file that can be used as a template
  # toolboxloc = location and name of folder containing cifti-matlab or fieldtrip toolbox
  # scratchloc = location where CSV files will be written then deleted
  
  # this function writes a MATLAB script and executes it
  # must have MATLAB installed
  
  # It 
  
  # Write table to CSV in scratch location
  fname <- file.path(scratchloc, 'tmp.csv')
  write.table(table, file=fname, row.names=FALSE, col.names=FALSE, sep=',')
  
  # Write MATLAB script
  line1 <- paste0("addpath '", toolboxloc1, "'")
  line2 <- paste0("addpath '", toolboxloc2, "'")
  line3 <- paste0("cifti = ciftiopen('",template,"','",wbloc,"');") 
  line4 <- paste0("data = csvread('",scratchloc,"/tmp.csv')",";")
  line5 <- "[nrow,ncol]=size(data);"
  
  # Zero out entires, allowing for files without subcortical voxels:
  line6 <- "cifti.cdata = NaN(91282,ncol);"
  line7 <- "cifti.cdata(1:nrow,:)=data;"
  line8 <- paste0("ciftisave(cifti,'",filename,"','",wbloc,"');")
  matlab_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8)
  writeLines(matlab_lines, con=file.path(scratchloc, 'myscript.m'))
  
  system(paste0(matlabloc," -nodisplay -r \"run('~/myscript.m'); exit\""))
  
  file.remove(fname)
  file.remove(file.path(scratchloc, 'myscript.m'))
}





##############
makeciftimmp <- function(v360, filename, template = '~/Dropbox/MyHCP/Data/100307/MNINonLinear/tfMRI_MOTOR_LR_Atlas.dtseries.nii', toolboxloc1 = '~/Dropbox/Applications2/hcp_ciftimatlab', toolboxloc2 = '~/Dropbox/mfunctions/robertoostenveld-cifti-matlab-27383b8',wbloc = '~/Applications2/workbench/bin_macosx64/wb_command',scratchloc = '~',matlabloc = '/Applications/MATLAB_R2018b.app/bin/matlab',mmpatlasloc="~/Dropbox/HCP_rsfmri_mmp/Data/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii"){
  
  # This function modifies makecifti to create a cifti file from 360 values corresponding to the glasser mmp atlas

    # Write table to CSV in scratch location
  fname <- file.path(scratchloc, 'tmp.csv')
  write.table(v360, file=fname, row.names=FALSE, col.names=FALSE, sep=',')
  
  # Write MATLAB script
  line1 <- paste0("addpath '", toolboxloc1, "'")
  line2 <- paste0("addpath '", toolboxloc2, "'")
  line3 <- paste0("cifti = ciftiopen('",template,"','",wbloc,"');") 
  line4 = paste0("mmpatlas = ciftiopen('",mmpatlasloc,"','",wbloc,"');")
  line5 = paste0("ncii = size(cifti.cdata,1);")
  line6 = paste0("ncortex = 59412;")
  line7 <- paste0("data = csvread('",scratchloc,"/tmp.csv')",";")
  line8 <- "[~,ncol]=size(data);"
  
  # Zero out entires, allowing for files without subcortical voxels:
  line9 <- "cifti.cdata = NaN(91282,ncol);"
  # TO DO: Expand to subcortical
  line10 <- "mmpatlas_plus = [mmpatlas.cdata;zeros(ncii-ncortex,1)];"
  line11 <- "for r = 1:360" 
  line12 <- "indices = mmpatlas_plus==r;"
  line13 <- "cifti.cdata(indices,:) = repmat(data(r,:),[sum(indices),1]);"
  line14 <- "end"
  line15 <- paste0("ciftisave(cifti,'",filename,"','",wbloc,"');")
  matlab_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15)
  writeLines(matlab_lines, con=file.path(scratchloc, 'myscript.m'))
  
  system(paste0(matlabloc," -nodisplay -r \"run('~/myscript.m'); exit\""))
  
  file.remove(fname)
  file.remove(file.path(scratchloc, 'myscript.m'))
}
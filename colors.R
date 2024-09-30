# cell type colors
cell_type_colors <- c(
    '#d8391f', '#bd6d84', '#1a5e62',
    '#81ab7d', '#213e87', '#71ad2c', 
    '#bc6a23', '#e2c4aa', '#6A3D9A',
    '#B15928', '#EA469A', "#FB9A99",
    '#88181d', '#7eb2d0', '#fcc02d',
    '#d8391f',
    '#6A3D9A', '#FB9A99', '#992CA0',
    '#213e87', '#FD9C6F', '#8ADF8D')

names <- c(
    'Epithelial_cancer', 'CAF', 'Macrophage',
    'Endothelial', 'T_cells', 'B_cells', 
    'DC', 'MAST', 'Myoepithelial',
    'PVL', 'Epithelial_basal', 'Epithelial_luminal',
    'Epithelial_KRT7+', 'Plasma', 'NK',
    "Epithelial_KRT7-",
    'Cholangiocyte', 'Hepatocyte', 'LSEC',
    'Mixed_lymphocytes', 'RBCs', 'Lymphatic_endothelial')


names(cell_type_colors) <- names
rm(names)


# old palette
# cell_type_colors <- c("#1F78B4", "#FFFF99", "#33A02C",
#     '#EA4648', '#6A3D9A', '#EA469A', '#FB9A99', '#6A3D9A',
#     '#B15928', '#FF7F00', '#A6CEE3', '#B2DF8A', '#CAB2D6', '#B15928',
#     '#FB9A99', '#992CA0', '#8ADF8D', '#FDBF6F', '#FD9C6F')
# names <- c('CAF', 'PVL', 'Endothelial', 
#     'Epithelial_cancer', 'Cholangiocyte', 'Epithelial_basal', 'Epithelial_luminal', 'Myoepithelial', 
#     'Mixed_lymphocytes', 'Macrophage', 'B_cells', 'DC', 'MAST', 'T_cells',
#     'Hepatocyte', 'LSEC', 'Lymphatic_endothelial', 'unknown', 'RBCs')




# snRNA-seq workflow colors
workflow_colors <- c("#0072B2", "#009E73", "#E69F00", '#00a3ff', '#CC79A7', '#a5682a')
workflow_names <- c('FFPE-snPATHO-Seq', 'Frozen-Flex', "Frozen-3'", 'FFPE-scFFPE', 'Gavish_et_al_2023', 'Spatial')
names(workflow_colors) <- workflow_names
rm(workflow_names)



# breast cancer sample colors
sample_colors <- c('#CC79A7', '#9e79cc', '#a7cc79', '#79cc9e')
sample_names <- c('Gavish_et_al_2023', '4411-mBC-Lum', '4066-pBC-HER2', '4399-mBC-TNBC')
names(sample_colors) <- sample_names
rm(sample_names)


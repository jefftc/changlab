/*clustering.pl*/
:- dynamic cluster_file/4.
/*file type
cluster_file(DatasetId,Contents,Parameters,Modules)
cluster_heatmap(DatasetId,Contents,Parameters,Modules)
class_neighbors(DatasetId,Contents,Parameters,Modules)
*/
%output interface
%cluster_file(Datasetid,Contents,Parameters,Modules):-
    
convert_cluster_parameters(Parameters,NewParameters):-
    convert_parameters_file(Parameters,NewParameters0),

    get_value_variable(Parameters,cluster_alg,Cluster_Alg),
    member(Cluster_Alg,[som,pca,kmeans,hierarchical]),
    append(NewParameters0,[cluster_alg,Cluster_Alg],NewParameters1),
    
    get_value_variable(Parameters,distance,Distance),
    member(Distance,[correlation,euclidean]),
    append(NewParameters1,[distance,Distance],NewParameters2),

    get_value(Parameters,k,5,K),
    not(atom(K)),
    append(NewParameters2,[k,K],NewParameters).

    
 /*-------------------------------------------------------------------------*/
% generate the cluster_file from signal_file with format pcl.    
cluster_file(DatasetId,Contents,Parameters,Modules):-
     convert_cluster_parameters(Parameters,NewParameters),
     get_value(NewParameters,format,unknown_format,Format),
     Format=pcl,
     convert_parameters_file(NewParameters,NewParameters1),
     signal_file(DatasetId,Contents,NewParameters1,Past_Modules),
     append(['DatasetId',DatasetId,'Contents',Contents],NewParameters,Write_list),
     Newadd = [clustering,Write_list],
     append(Past_Modules,Newadd,Modules).

 /*-------------------------------------------------------------------------*/

% generate the cluster_heatmap from cluster_file
cluster_heatmap(DatasetId,Contents,Parameters,Modules):-
    get_value_variable(Parameters,color,Color),
    member(Color,[blue_yellow,red_green]),
    (convert_cluster_parameters(Parameters,NewParameters),
    cluster_file(DatasetId,Contents,NewParameters,Past_Modules),
    append(NewParameters,[color,Color],NewParameters1),
    append(['DatasetId',DatasetId,'Contents',Contents],NewParameters1,Write_list);
    get_value_variable(Parameters,cluster_alg,Cluster_Alg),
    Cluster_Alg=no_cluster_alg,
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    signal_file(DatasetId,Contents,NewParameters,Past_Modules),
    append(NewParameters,[color,Color],NewParameters1),
    append(['DatasetId',DatasetId,'Contents',Contents,cluster_alg,
          no_cluster_alg],NewParameters1,Write_list)),
    Newadd = [cluster_heatmap,Write_list],
    append(Past_Modules,Newadd,Modules).

/*-------------------------------------------------------------------------*/
class_neighbors(DatasetId,Contents,Parameters,Modules):- 
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=gct,
    signal_file(DatasetId,Contents,NewParameters,Past_Modules_1),
    class_label_file(DatasetId,Contents,_,Past_Modules_2),
    append(Past_Modules_1,Past_Modules_2,Past_Modules),
    append(['DatasetId',DatasetId,'Contents',Contents],NewParameters,Write_list),
    Newadd = [class_neighbors,Write_list],
    append(Past_Modules,Newadd,Modules).

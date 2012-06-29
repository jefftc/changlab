/*clustering.pl*/
:- dynamic cluster_file/2.
/*file type
cluster_file(Parameters,Modules)
cluster_heatmap(Parameters,Modules)
class_neighbors(Parameters,Modules)
*/
%output interface

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
cluster_file(Parameters,Modules):-
     convert_cluster_parameters(Parameters,NewParameters),
     get_value(NewParameters,format,unknown_format,Format),
     Format=pcl,
     get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
     Is_Logged=logged,
     convert_parameters_file(NewParameters,NewParameters1),
     signal_file(NewParameters1,Past_Modules),
     Newadd = [clustering,NewParameters],
     append(Past_Modules,Newadd,Modules).

 /*-------------------------------------------------------------------------*/

% generate the cluster_heatmap from cluster_file
cluster_heatmap(Parameters,Modules):-
    (convert_cluster_parameters(Parameters,NewParameters1),
    cluster_file(NewParameters1,Past_Modules);
    get_value_variable(Parameters,cluster_alg,Cluster_Alg),
    Cluster_Alg=no_cluster_alg,
    convert_parameters_file(Parameters,NewParameters),
    get_value(NewParameters,format,unknown_format,Format),
    Format=pcl,
    get_value(NewParameters,is_logged,unknown_logged,Is_Logged),
    Is_Logged=logged,
    signal_file(NewParameters,Past_Modules),
    append(NewParameters,[cluster_alg,
          no_cluster_alg],NewParameters1)),
    get_options(Parameters,[hm_width,hm_height,color],[],Options),
    append(NewParameters1,Options,NewParameters2),
    Newadd = [make_heatmap,NewParameters2],
    append(Past_Modules,Newadd,Modules).



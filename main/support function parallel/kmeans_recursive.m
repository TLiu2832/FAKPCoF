function [Centre, minCentre]=kmeans_recursive(X,Cluster)
% function Centre=kmeans_recursive(X,Cluster)
%   kmeans_recursive    Bisecting k-means clustring
%   [IDX, C] = kmeans_recursive(X, delta) partititions the N x P data matrix X into
%   clusters through a fully vectorized algorithm, where N is the number of data points and 
%   P is the number of dimensions (variables). First whole data is bisected into clusters .
%   Chooses the cluster with point-centre distance is maximum and bisect it again to two clusters
%   This is iterated till point-centre distance is less than delta
%   IDX is the returned N x 1 vector contains the cluster indices of each point.
%   C is the K cluster centroid locations in the K x P matrix.
%
% Initializing cluster count
K=1;                                                                
% Performing K means clustering
[minCentre,Centretemp,newdistvect]=kmeans_cluster(X);       

d=size(Centretemp,2);
Centre=zeros(Cluster,d);
var=zeros(Cluster,1);

K=K+1;
% Incrementing cluster count 
Centre(1,:)=Centretemp(1,:); 
Centre(K,:)=Centretemp(2,:);    

var(1)=newdistvect(1); 
var(K)=newdistvect(2);
% Calculating maximum point-centre distance and cluster to be bisected  
[~,maxindex]=max(var);                                    
while (K<Cluster)    
    % Performing K means clustering   
    [label,Centretemp,newdistvect,newclustvect]=kmeans_cluster(X(minCentre==maxindex,:));     
   if (newclustvect(1) == 0)
   var(maxindex)=0; 
   elseif (newclustvect(2) == 0)
   var(maxindex)=0; 
   else
    % Incrementing cluster count
    K=K+1;                                                       
        
    % Populating new cluster centres
    Centre(maxindex,:)=Centretemp(1,:); 
    Centre(K,:)=Centretemp(2,:);   
  
    % Populating new cluster indices vector
    label=(label==1).*maxindex+(label==2).*K;
    minCentre(minCentre==maxindex)=label;   
        
    var(maxindex)=newdistvect(1); 
    var(K)=newdistvect(2);
    
   end
   if(any(var)==0) 
       break; 
   end
    % Calculating maximum point-centre distance and cluster to be bisected
    [~,maxindex]=max(var);                                
end
end

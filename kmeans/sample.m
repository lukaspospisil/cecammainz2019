clc
clear all

%close all
addpath('ProgramFiles/')
addpath('Plotting/')
addpath('data/')

% set random generator
randn('seed',13);
rand('seed',13);

K=3;

% problem parameters
T = 1e4; % number of points
N_anneal = 10; % number of annealing steps
myeps = 1e-4; % precision

X = generate_problem(T,2);

if false
    % set to true if you don't want to test GPU
    % don't forget to compile "CudaFiles/mykernel.cu" using "nvcc mykernel.cu -ptx"
    gpu = gpuDevice(1);
    reset(gpu);
    
    [out_kmeans] = mykmeans_gpu(X, K, N_anneal,myeps);
    wait(gpu)
    
else
    [out_kmeans] = mykmeans_cpu(X, K, N_anneal,myeps);
end

plot_classification( X,out_kmeans.gamma, 1:2)


N = 100;  %number of neurons
num_of_iterations = 2000;
threshold = 400;

%creating a vector from image
Img1 = imread("C:\Users\or\Documents\6.png");
newImg1 = imresize(Img1,[10,10]);
Img2 = imread("C:\Users\or\Documents\3.png");
newImg2 = imresize(Img2,[10,10]);
Img3 = imread("C:\Users\or\Documents\9.png");
newImg3 = imresize(Img2,[10,10]);

surf(peaks);

colormap(gray);



% creating 3 memories randomly
mem_1 = randi([0,1],10);
mem_2 = randi([0,1],10);
mem_3 = randi([0,1],10);

for c = 1:10
    for r = 1:10
        if newImg1(c,r)>100
            mem_1(c,r) = 1;
        else
            mem_1(c,r) = 0;
        end
    end
end
for c = 1:10
    for r = 1:10
        if newImg3(c,r)>100
            mem_3(c,r) = 1;
        else
            mem_3(c,r) = 0;
        end
    end
end
for c = 1:10
    for r = 1:10
        if newImg3(c,r)>100
            mem_3(c,r) = 1;
        else
            mem_3(c,r) = 0;
        end
    end
end

% setting memories to be vectors
mem1_vec = reshape(mem_1,[],1);
mem2_vec = reshape(mem_2,[],1);
mem3_vec = reshape(mem_3,[],1);

%combinig to a matrix of 3X100
memories = [mem1_vec,mem2_vec,mem3_vec];


%showing that for different thresholds we have different success percent
thresholds = [400, 400];
number_of_succ = zeros(1,5);
converge_ave =zeros(1,5);
for i=1:1
    successes =0;
    conver = 0;
    for j = 1:1
    [suc,conv,per_of_simil] = one_run(memories,N,num_of_iterations,thresholds(i));   %run function of one attempt of Hopfield model
    disp(['converged in step: ',num2str(conv)])
    conver = conver + conv/100;
    if suc == 1
        disp('to the right memory')
        successes = successes +1; 
    else
        disp('not to the right memory')
    end 
    end
    converge_ave(i) = conver;
    number_of_succ(i) = successes;
end

figure
scatter(thresholds,number_of_succ/100)

%showing that for each run we have different convergence time

% try_num = 1:20;
% converence = zeros(1,20);
% for j = 1:20
%     [suc,conv,per_of_simil] = one_run(memories,N,num_of_iterations,threshold);   %run function of one attempt of Hopfield model
%     disp(['converged in step: ',num2str(conv)])
%     converence(j) = conv;
% end
% 
% figure
% bar(try_num,converence)

% %show plot for 1 convergence
% [suc,conv,per_of_simil] = one_run(memories,N,num_of_iterations,threshold);   %run function of one attempt of Hopfield model
% figure
% plot((1:conv),per_of_simil(1:conv))


function [success,converge,percent_of_similarity] = one_run(memories,N,number_of_iter,threshold)

A = set_matrix_of_memories(memories);
w = create_weights(N,A);
noise = 0.1;  %noise of 10 percent
rand_mem_ind = 3;  %randomly choosing a memory to change
rand_mem = A(:,rand_mem_ind); 
noise_ind = randperm(length(rand_mem),noise*length(rand_mem));  %randomly choosing bits indexes to change
noise_mem = rand_mem;
noise_mem(noise_ind) = rand_mem(noise_ind)*(-1);  %making noise
m1 = noise_mem;
mat = reshape(m1, 10,10);
for c = 1:10
    for r = 1:10
        if mat(c,r)>0
            mat(c,r) = 1000;
        else
            mat(c,r) = 0;
        end
    end
end
saveas(image(mat),'z.png');
pause(7);
s_0 = noise_mem;  %defining the maximum bumber of iterations for the system
s=zeros(N,number_of_iter);  %defining a matrix containing the dynamics of the system
% meaning each column i is the state of the system in iteration i
s(:,1) = s_0;   %setting the initial state to be the noise_memory


count_no_change = 0;
count_diff = 0;
converge = 0;
success = 0;

percent_of_similarity = zeros(1,number_of_iter);

for i=1:number_of_iter 
    rand_neuron = randi([1,100],1);  %in order to work asynchronous we choose one neuron randomly 
    h_i = 0;
    for j=1:N  %calculate all inputs from other neurons
        h_i = h_i + w(rand_neuron,j)*s(j,i);
    end
    s_i1 = sign(h_i);  %calculate the state of the random neuron in the next step
    s(:,i+1) = s(:,i);
    s(rand_neuron,i+1) = s_i1;  %setting the whole's system state for the next step
    if mod(i,10) == 0
    tmp = reshape(s(:,i), 10,10);
        for c = 1:10
            for r = 1:10
                if tmp(c,r)>0
                    tmp(c,r) = 1000;
                else
                    tmp(c,r) = 0;
                end
            end
        end
        saveas(image(tmp),'z.png')
    end

% % for convergence graph: this is not part of the algorithm, it is only to
% determine how close the system is to the real solution at that iteration
    how_similiar = 0;
    for l=1:N
        if s(l,i+1) == rand_mem(l)
            how_similiar = how_similiar + 1;
        end
    end
    percent_of_similarity(1,i) = how_similiar/100;
    
    
    %here i create a counter that countes how many steps the system has
    %been stable
    if s(:,i+1) == s(:,i)
        count_no_change = count_no_change + 1;
%         disp 'same'
    else
%         disp 'diff'
        count_no_change = 0;
        count_diff = count_diff + 1;
    end
    
    %if the number of steps the system has been stable is above certain
    %value we determine convergenes, and also check if it converged to the
    %right memory
    if count_no_change >= threshold
        converge = i;
%         disp (['converged for sure in iteration: ', num2str(i)])
%         disp (['there are: ' , num2str(count_diff), ' different ones'])
        converged_mem = s(:,converge);
        if converged_mem == rand_mem
%             disp 'and to the right memory!'
            success = 1;
            final = reshape(s(:,i), 10,10);
            for c = 1:10
                for r = 1:10
                    if final(c,r)>0
                        final(c,r) = 1000;
                    else
                        final(c,r) = 0;
                    end
                end
            end
            saveas(image(final),'z.png')
        
            
        else
%             disp 'but not to the right memory'
        end
        break
    end 
end
end


function w = create_weights(N,A)
    w = (1/N)*(A*A');
    for k=1:N
    for l=1:N
        if k==l
            w(k,l) = 0;
        end
    end
    end
end

function A = set_matrix_of_memories(memories)
A = memories;
for i=1:100
    for j=1:3
        if A(i,j) == 0
            A(i,j) = -1;
        end
    end
end
end

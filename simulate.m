function [ bcost_mon, bcost_ntfy, bcost_total, dcost_mon, dcost_ntfy, dcost_qry, dcost_total ] = simulate( tau, duration, lmd, p, c, cntfy, f )

%------- MATLAB INPUT ----------
%tau = 1/1;
%duration = 10000;
%lmd = [ 0.0117170204491856 0.0147226586946997 0.0114397637358739 0.0108102453512919 0.0139518577451191 0.0152237583656539 ];
%p = [ 0.308300779743761 0.389918658212339 0.335335925669604 0.504876905354311 0.58870212697272 0.517017747204504 ];
%c = [ 0.599521553015445 0.604212284879583 0.681357338477975 0.545010426731755 0.667109649005177 0.449336304020588 ];
%cntfy = [ 0 0 0 0 0 0 ];
%F = [ 
%1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
%2 1 2 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
%2 2 4 5 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
%];
%cntfy = [ 0 0 0 0 0 ] % as many as subfs


% simulation
%simSmartMonitor()

%cost of basic

format long g

% fix: use composite lambda
sn=length(lmd);
bcost_mon = duration * tau * sum(c);
%bcost_ntfy = duration * cntfy * lambda_sf';
g=sprintf('%d ',lmd);
fprintf('LMD: %s\n', g);
g=sprintf('%d ',p);
fprintf('P: %s\n', g);
g=sprintf('%d ',c);
fprintf('C: %s\n', g);



[m, n] = size(f);
indexi = 1;
indexj = 1;
F = zeros((n/31), 31);
for i=1:n
    if mod(i,31) == 0
        F(indexi, indexj) = f(i);
        indexi = indexi+1;
        indexj = 1;
        continue;
    end
    F(indexi, indexj) = f(i);
    indexj = indexj + 1;
end

%disp(F);
% TODO
% notifying cost of basic
% notify whenver sf changes
% bcost_ntfy = sum_{f:subf} \lambda_f
% \lambda_f = 
bcost_ntfy = 0;
dcost_ntfy = 0;


mcost_cl_lst=[]; qcost_cl_lst=[]; qcc_cl_lst=[];
qcost_total = 0; prb_cl_lst = []; lambda_cl_lst = [];
%for each clause
for i=1:size(F,1)
    nf = F(i,1);
    cur=2; % current position in F's row
    
    atoms = [];
    % for each subformula
    for j = 1: F(i,1)
        na = F(i,cur); % number of atoms in this subformula
        atoms = [atoms; F(i,cur+1:cur+na)'];
        cur = cur + na + 1;
    end
    %g=sprintf('%d ',atoms);
    %fprintf('Atoms: %s\n',g);
    psubf = p(atoms);
    lmdsubf = lmd(atoms);
    csubf = tau*c(atoms);
    %g=sprintf('%d ',psubf);
    %fprintf('PSUBF: %s\n',g);
    %g=sprintf('%d ',lmdsubf);
    %fprintf('LMDSUBF: %s\n',g);
    %g=sprintf('%d ',csubf);
    %fprintf('CSUBF: %s\n',g);
    
    % cost of monitoring i'th clause and its querying overhead
    [mcost_cl, ppi_cl, qcost_cl, qcc_cl] = costdim( csubf', psubf', lmdsubf', csubf' ); 
    CL_pi = ppi_cl;
    CL_monprob = ppi_cl(1:length(ppi_cl)-1) + ppi_cl(length(ppi_cl));
    prb_cl_lst = [ prb_cl_lst; prod(psubf) ];
    s=0;
    len = length(lmdsubf);
    for k=1:len
        s = s + prod(psubf(1:k-1))*prod(psubf(k+1:len))*lmdsubf(k);
    end
    lambda_cl_lst= [lambda_cl_lst; s]; % change rate of the clause
    mcost_cl_lst = [ mcost_cl_lst; mcost_cl ];
    qcost_cl_lst = [qcost_cl_lst; qcost_cl];
    qcc_cl_lst = [ qcc_cl_lst; qcc_cl ];
%    basic_cl = [basic_cl; sum(basic_f) ];
end
qcost_total = qcost_total + sum( qcost_cl_lst );

[mcost_be, ppi_be, qcost_be, qcc_be] = costdim( mcost_cl_lst, 1-prb_cl_lst, lambda_cl_lst, qcc_cl_lst );
ppi_be;
bcost_total = bcost_mon + bcost_ntfy;

dcost_mon=mcost_be*duration;
dcost_qry=qcost_be*duration;
dcost_total = dcost_mon + dcost_qry;

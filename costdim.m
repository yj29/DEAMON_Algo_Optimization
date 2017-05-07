function [costresult, pie, costq, qcc] = costdim(carr,parr,larr, qccarr)
%carr; costs, parr: probs, larr: lambdas
n=length(carr); %number of elements
o=n+1; %refers state zero

if n == 1
    costresult = carr(1);
    pie = [ 1; 0 ];
    costq = 0;
    qcc = carr(1);
    return
end

lmd1 = larr./(1-parr)./2;
lmd0 = larr./parr./2;

%compute transition matrix
%for each state 1,2,...,n+1
for i=1:n
%    lmd1 = larr(i)/(2*(1-parr(i)));
%    lmd0 = larr(i)/(2*parr(i));
    %lmd1 = larr(i)*parr(i);
    %lmd0 = larr(i)*(1-parr(i)); 
    for j=1:n
        if i == j
            q(i,j) = 0;
            continue;
        end
        q(i,j) = lmd1(i) * (1-parr(j));
        for k=1:j-1
            if k ~= i
                q(i,j) = q(i,j) * parr(k);
            end
        end
        %q(i,j)=lmd1*prod(parr(1:j-1))*(1-parr(j));
    end
    q(i,o)=lmd1(i);
    for k=1:n
        if k ~= i
            q(i,o) = q(i,o) * parr(k);
        end
    end
    q(o,i)=lmd0(i)/2;
   % q(o,i)=sum(lmd0)/length(lmd0);
    %q(o,i) = sum(1./lmd0)/length(lmd0)
end
for i=1:o
    q(i,i)=-sum(q(i,:));
end
qp=[q ones(o,1)];
b=[zeros(o,1);1];
q;
pie=qp'\b;
costresult=pie'*[carr; sum(carr)];

% now compute querying cost: one time query to determine next state

% build query cost matrix
% qcost = sum_{i,j (i!= j) pie_i q_{i,j} qc_{i,j}
% qc_{i,j} = sum_{k=1,k!=j}^{j} c_k
cpre = tril(ones(n+1,n+1))*[qccarr; 0];
qcc = 0; % average cost of querying this when disabled
for i=1:n+1
    for j=1:n+1
        if i == j
            qc(i,j) = 0;
        elseif i == n+1
            qc(i,j) = 0;
        else
            qc(i,j) = cpre(j);
            if (i < n+1) && (i < j)
                qc(i,j) = qc(i,j)-qccarr(i);
            end
            qc(i,j) = qc(i,j) * pie(i);
        end
    end
    
    if i == n+1
        1;
    else
        qcc = qcc + prod(parr(1:i-1)) * (1-parr(i)) * cpre(i);
    end
end

% sum of costs for querying elements during monitoring
costq = sum(sum(q.*qc));




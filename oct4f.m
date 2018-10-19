
function dy = oct4f(~,y,logp)

p = 10.^logp;

% species 
diff_signal=y(1);
oct4_mrna=y(2);
oct4_prot=y(3);
cdx2_mrna=y(4);
cdx2_prot=y(5);

% parameters
K1 = p(1);
K2 = p(2);
n1 = p(3);
n2 = p(4);
beta1 = p(5);
beta2 = p(6); 

%degradation constant.
a1 = p(7);
% a1 = log(2)/10;
a2 = log(2)/8.65;
a3 = p(8);
a4 = p(9);
s1 = p(10);
s2 = p(11);
% The following are the numerator and denominator terms for OCT4 mrna fc
oct4fcnum = a1/(1+diff_signal+(cdx2_prot/K1)^n1);
oct4fcden=1/(1+(cdx2_prot/K1)^n1);
% oct4 steady state
oct4_st = (s1/a1)*(1/1+diff_signal+(1/K1)^n1);
% OCT4 protein over steady state gives us fc of OCT4 protein
oct4mRNA_prod = oct4fcnum/oct4fcden;
% OCT4 mrna degradation
oct4mRNA_deg = (a1*oct4_mrna)/oct4_st;
% OCT4 protein production and degradation 
oct4prot_prod = beta1*oct4_mrna;
oct4prot_deg = a2*oct4_prot;

% CDX2 mrna fold change numerator and denominator
cdx2num = a3/(1+(oct4_prot/K2)^n2);
cdx2dem=1/(1+(oct4_prot/K2)^n2);
cdx2mRNA_prod = cdx2num/cdx2dem;
% cdx2 steady state
cdx2_st = (s2/a3)*(1/1+(1/K2)^n2);

% CDX2 degradation
cdx2mRNA_deg = a3*cdx2_mrna/cdx2_st;

% CDX2 protein production and degradation
cdx2prot_prod = beta2*cdx2_mrna;
cdx2prot_deg = a4*cdx2_prot;

% production and degradation of the differentiation signal
diff_prod = 0;
diff_deg = 0;

% fold change difference of OCT4 and CDX2
dy = [diff_prod - diff_deg;
      oct4mRNA_prod-oct4mRNA_deg;
      oct4prot_prod - oct4prot_deg;
      cdx2mRNA_prod-cdx2mRNA_deg;
      cdx2prot_prod - cdx2prot_deg;
     ];
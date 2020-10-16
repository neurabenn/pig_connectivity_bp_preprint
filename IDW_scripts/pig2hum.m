
function hmap  = pig2hum(hum_BP,pig_BP,human,pig,varargin)
%%%% inouts are a human BP, a pig BP of equal size, a pig .func.gii and the
%%%% human surf.gii 
[~,BP_out,~]=fileparts(hum_BP);
pig_BP=load(pig_BP);
pig_BP=pig_BP.bp;
hum_BP=load(hum_BP);
hum_BP=hum_BP.bp; %hum_BP_29;
pig_srf=gifti(pig);
hum_srf=gifti(human);


pig_map=pig_srf.cdata;
pig_map(pig_map==0)=0.0000000000001;

%%%interpolate
N=10001;
gamma = -4;
PH=calc_KL(hum_BP,pig_BP);
% interpolation - careful with division by zero etc.
n     = size(PH,2);
D     = PH .* repmat(~~pig_map',n,1);
D     = D.^gamma;  D(isnan(D))=0; D(isinf(D))=0;
W     = D ./ repmat(sum(D,2)+~sum(D,2),1,n);
%%%% do interpolatoin as matrix multiplication. Save gifti out
hmap=W*pig_map;

hmax=max(hmap);
if length(varargin)>0
    thr=varargin{1}*hmax;
    hmap(hmap<thr)=0;
end

hmap=gifti(hmap);
%%%% set out put path
[filepath,name] = fileparts(pig);
% filepath='/Volumes/SC/LAB_BI/LAB/Austin_Benn/pig_BP_paper/figures/surfICAS/left_surf';
out=sprintf('%s/%s_%s.gii',filepath,BP_out,name);
display(out)
save(hmap,out,'Base64Binary');

end
#### adapted from RMARS https://doi.org/10.7554/eLife.35237


function pmap  = hum2pig(hum_BP,pig_BP,human,pig,varargin)
%%%% inouts are a human BP, a pig BP of equal size, a pig.surf.gii and the
%%%% human func.gii 
[~,BP_out,~]=fileparts(hum_BP);
pig_BP=load(pig_BP);
pig_BP=pig_BP.bp;
hum_BP=load(hum_BP);
hum_BP=hum_BP.bp; %hum_BP_29;
pig_srf=gifti(pig);
hum_srf=gifti(human);


hum_map=hum_srf.cdata;
hum_map(hum_map==0)=0.0000000000001;

%%%interpolate
N=10001;
gamma = -4;
PH=calc_KL(pig_BP,hum_BP);
% interpolation - careful with division by zero etc.
n     = size(PH,2);
D     = PH .* repmat(~~hum_map',n,1);
D     = D.^gamma;  D(isnan(D))=0; D(isinf(D))=0;
W     = D ./ repmat(sum(D,2)+~sum(D,2),1,n);
%%%% do interpolatoin as matrix multiplication. Save gifti out
pmap=W*hum_map;

pmax=max(pmap);
if length(varargin)>0
    thr=varargin{1}*hmax;
    pmap(hmap<thr)=0;
end

pmap=gifti(pmap);
%%%% set out put path
[filepath,name] = fileparts(human);
filepath='/Volumes/SC/LAB_BI/LAB/Austin_Benn/pig_BP_paper/figures/surfICAS/left_surf';
out=sprintf('%s/%s_2pig_%s.gii',filepath,BP_out,name);
save(pmap,out,'Base64Binary');

end
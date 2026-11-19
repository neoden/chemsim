const fs=require('fs');
const {build,solve}=new Function(fs.readFileSync(__dirname+'/eht-core.js','utf8')+';return {build,solve};')();
function chi(a,t,p){const r=[p[0]-a.r[0],p[1]-a.r[1],p[2]-a.r[2]],r2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];let s=0;const c=t?a.dp:a.ds;for(let i=0;i<3;i++)s+=c[i]*Math.exp(-a.al[i]*r2);return t?s*r[t-1]:s;}
for(const name of ['h2o','nh3','bh3','c2h4','pair']){
  const atoms=build(name,1.66),R=solve(atoms);
  console.log('\n==',name,'basis',R.n,'nocc',R.nocc);
  console.log(' E:',R.E.map(e=>e.toFixed(3)).join(' '));
  // residual HC-SCE and C^T S C
  let res=0,orth=0;
  for(let m=0;m<R.n;m++){for(let i=0;i<R.n;i++){let hc=0,sc=0;for(let j=0;j<R.n;j++){hc+=R.H[i][j]*R.C[m][j];sc+=R.S[i][j]*R.C[m][j];}res=Math.max(res,Math.abs(hc-R.E[m]*sc));}
    for(let m2=0;m2<R.n;m2++){let t=0;for(let i=0;i<R.n;i++)for(let j=0;j<R.n;j++)t+=R.C[m][i]*R.S[i][j]*R.C[m2][j];orth=Math.max(orth,Math.abs(t-(m===m2?1:0)));}}
  console.log(' max|HC-SCE|',res.toExponential(2),' max|CtSC-I|',orth.toExponential(2));
  const desc=m=>R.C[m].map((c,i)=>Math.abs(c)>0.15?`${atoms[R.bs[i].k].el}${R.bs[i].k}${'spxpypz'.match(/s|p./g)[R.bs[i].t]}:${c.toFixed(2)}`:'').filter(Boolean).join(' ');
  console.log(' HOMO',R.E[R.nocc-1].toFixed(3),desc(R.nocc-1));
  console.log(' LUMO',R.E[R.nocc].toFixed(3),desc(R.nocc));
}
// numeric integration check on water: analytic S vs grid, and total electrons
{
  const atoms=build('h2o'),R=solve(atoms),n=R.n,h=0.1,L=8,Sn=[];for(let i=0;i<n;i++)Sn.push(new Array(n).fill(0));
  let ne=0;const v=new Array(n);
  for(let x=-L;x<=L;x+=h)for(let y=-L;y<=L;y+=h)for(let z=-L;z<=L;z+=h){
    for(let i=0;i<n;i++)v[i]=chi(atoms[R.bs[i].k],R.bs[i].t,[x,y,z]);
    for(let i=0;i<n;i++)for(let j=i;j<n;j++)Sn[i][j]+=v[i]*v[j];
    for(let m=0;m<R.nocc;m++){let q=0;for(let i=0;i<n;i++)q+=R.C[m][i]*v[i];ne+=2*q*q;}
  }
  let dmax=0;for(let i=0;i<n;i++)for(let j=i;j<n;j++)dmax=Math.max(dmax,Math.abs(Sn[i][j]*h*h*h-R.S[i][j]));
  console.log('\nwater grid check: max|S_num-S_analytic| =',dmax.toExponential(2),' electrons =',(ne*h*h*h).toFixed(4));
}

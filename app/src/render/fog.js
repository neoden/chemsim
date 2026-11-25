// Volumetric "fog" renderer: WebGL2 ray marching in the fragment shader.
// The basis set and the orbital coefficients come in as data textures, so the number of
// atoms, shells and primitives per shell is whatever the engine reports.
//
// Looks:
// absorption along the ray, soft translucent shells, lighting by the analytic gradient,
// two colours by the sign of the wavefunction, dim grey "body" of the total valence density,
// pale "ghost" for empty orbitals, nuclei as small glowing dots.

const VS = `#version 300 es
void main(){vec2 p=vec2(float((gl_VertexID<<1)&2),float(gl_VertexID&2));gl_Position=vec4(p*2.0-1.0,0.0,1.0);}`;

// Every fragment shader here declares `precision highp int`. Ints are mediump by default
// there, and mediump float->int conversions (texel -> index) crash Mesa's ACO shader compiler
// on some Radeon cards (radeonsi: "Unknown NIR ALU instr: f2fmp"), which kills the browser's
// GPU process and loses the WebGL context.

// Shared GLSL: the basis set as data textures.
const BASIS_GLSL = `
// uAtom: row 0 = xyz, cutoff radius^2; row 1 = first shell, shell count, fragment (0 or 1)
// uShell: l, first primitive, primitive count, first AO. A shell has 2l+1 functions, l <= 2:
// s; p: x, y, z; d: sqrt(3) xy, sqrt(3) yz, z2-(x2+y2)/2, sqrt(3) xz, sqrt(3)/2 (x2-y2)
// uPrim: exponent, coefficient
// uCoef: the selected orbital
// uOcc: the occupied orbitals times sqrt(occupation), four to a row (one in each channel),
// uNO rows: the density sums four orbitals with every texel it reads
uniform sampler2D uAtom,uShell,uPrim,uCoef,uOcc;
uniform int uNA,uNO;
const float S32=0.8660254038; // sqrt(3)/2
float cf(int ao,int row){return texelFetch(uCoef,ivec2(ao,row),0).x;}
// selected orbital with its analytic gradient, split into the parts that sit on the atoms
// of fragment 0 (pa, ga) and of fragment 1 (pb, gb); psi = pa + pb
void orbital2(vec3 p,out float pa,out vec3 ga,out float pb,out vec3 gb){
pa=0.0;pb=0.0;ga=vec3(0.0);gb=vec3(0.0);
for(int a=0;a<uNA;a++){
vec4 A=texelFetch(uAtom,ivec2(a,0),0);
vec3 r=p-A.xyz;
float r2=dot(r,r);
if(r2>A.w)continue; // far atom: its functions are negligible here
vec3 I3=texelFetch(uAtom,ivec2(a,1),0).xyz;
ivec2 I=ivec2(I3.xy);
float psi=0.0;vec3 g=vec3(0.0);
for(int s=I.x;s<I.x+I.y;s++){
ivec4 S=ivec4(texelFetch(uShell,ivec2(s,0),0));
float R=0.0,R1=0.0;
for(int k=S.y;k<S.y+S.z;k++){
vec2 q=texelFetch(uPrim,ivec2(k,0),0).xy;
float e=q.y*exp(-q.x*r2);
R+=e;R1-=2.0*q.x*e;
}
if(S.x==0){
float c=cf(S.w,0);
psi+=c*R;g+=c*R1*r;
}else if(S.x==1){
vec3 c=vec3(cf(S.w,0),cf(S.w+1,0),cf(S.w+2,0));
float cr=dot(c,r);
psi+=cr*R;g+=cr*R1*r+c*R;
}else{
// d: the five functions together are a quadratic form, r.Qr times R; q = Qr
float c0=cf(S.w,0),c1=cf(S.w+1,0),c2=cf(S.w+2,0),c3=cf(S.w+3,0),c4=cf(S.w+4,0);
vec3 q=vec3(-0.5*c2+S32*c4,-0.5*c2-S32*c4,c2)*r+S32*vec3(c0*r.y+c3*r.z,c0*r.x+c1*r.z,c3*r.x+c1*r.y);
float qr=dot(q,r);
psi+=qr*R;g+=qr*R1*r+2.0*q*R;
}
}
if(I3.z<0.5){pa+=psi;ga+=g;}else{pb+=psi;gb+=g;}
}
}
void orbital(vec3 p,out float psi,out vec3 g){
float pa,pb;vec3 ga,gb;
orbital2(p,pa,ga,pb,gb);
psi=pa+pb;g=ga+gb;
}
// total valence density: sum over occupied orbitals
float density(vec3 p){
float chi[MAXAO];int idx[MAXAO];int n=0;
for(int a=0;a<uNA;a++){
vec4 A=texelFetch(uAtom,ivec2(a,0),0);
vec3 r=p-A.xyz;
float r2=dot(r,r);
if(r2>A.w)continue;
ivec2 I=ivec2(texelFetch(uAtom,ivec2(a,1),0).xy);
for(int s=I.x;s<I.x+I.y;s++){
ivec4 S=ivec4(texelFetch(uShell,ivec2(s,0),0));
float R=0.0;
for(int k=S.y;k<S.y+S.z;k++){
vec2 q=texelFetch(uPrim,ivec2(k,0),0).xy;
R+=q.y*exp(-q.x*r2);
}
if(S.x==0){chi[n]=R;idx[n]=S.w;n++;}
else if(S.x==1){chi[n]=r.x*R;idx[n]=S.w;chi[n+1]=r.y*R;idx[n+1]=S.w+1;chi[n+2]=r.z*R;idx[n+2]=S.w+2;n+=3;}
else{
float t=2.0*S32*R;
chi[n]=t*r.x*r.y;chi[n+1]=t*r.y*r.z;chi[n+2]=(r.z*r.z-0.5*(r.x*r.x+r.y*r.y))*R;chi[n+3]=t*r.x*r.z;chi[n+4]=S32*R*(r.x*r.x-r.y*r.y);
for(int j=0;j<5;j++)idx[n+j]=S.w+j;
n+=5;
}
}
}
float rho=0.0;
for(int m=0;m<uNO;m++){
vec4 q=vec4(0.0);
for(int j=0;j<n;j++)q+=texelFetch(uOcc,ivec2(idx[j],m),0)*chi[j];
rho+=dot(q,q);
}
return rho;
}`;

// The grids are filled in one draw call each: all the slices of the cube side by side on a
// flat sheet, uCols of them in a row, which is then copied into the 3D texture slice by slice.
// One draw call a slice was several times slower for the body: a call takes at least as long
// as the shader runs for one point, however few points there are in it.
const SHEET_GLSL = `
uniform float uRb,uN,uCols;
// the point of the cube that this fragment of the sheet stands for; false past the last slice
bool sheetPoint(out vec3 p){
vec2 tile=floor(gl_FragCoord.xy/uN);
float z=tile.y*uCols+tile.x;
p=(vec3(gl_FragCoord.xy-tile*uN,z+0.5)/uN*2.0-1.0)*uRb;
return z<uN-0.5;
}`;

// Fills the sheet with log(rho) for the body-density grid.
const gridSource = (maxAo) => `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
#define MAXAO ${maxAo}
out vec4 o;
${SHEET_GLSL}
${BASIS_GLSL}
void main(){
vec3 p;
if(!sheetPoint(p))discard;
o=vec4(log(density(p)+1e-13),0.0,0.0,1.0);
}`;

// Fills the sheet for the orbital grid: the selected orbital with its analytic gradient, the
// part that sits on the atoms of fragment uPart.
const orbitalGridSource = (maxAo) => `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
#define MAXAO ${maxAo}
out vec4 o;
uniform int uPart;
${SHEET_GLSL}
${BASIS_GLSL}
void main(){
vec3 p;
if(!sheetPoint(p))discard;
float pa,pb;vec3 ga,gb;
orbital2(p,pa,ga,pb,gb);
o=uPart==0?vec4(pa,ga):vec4(pb,gb);
}`;

// Colour of the cut at point p: the field is evaluated exactly there (no grid), drawn as a
// colour map with isolines every factor of 2 and, for an orbital, the nodal line in white.
const CUT_GLSL = `
const vec3 CPc=vec3(0.22,0.52,1.0);
const vec3 CNc=vec3(1.0,0.5,0.12);
uniform int uPlField; // 0: selected orbital psi, 1: total density (shown as sqrt(rho))
float cutField(vec3 p){
if(uPlField==0){float psi;vec3 g;orbital(p,psi,g);return psi;}
return sqrt(density(p));
}
vec4 cutColor(vec3 p){
float v,gate=0.0;
if(uPlField==0){
float psi;vec3 g;orbital(p,psi,g);v=psi;
gate=smoothstep(1.5e-3,4e-3,length(g)); // no nodal line far out, where psi is just numerically zero
}else v=sqrt(density(p));
vec3 base=uPlField==0?(v>0.0?CPc:CNc):vec3(0.74,0.77,0.84);
float m=abs(v)+1e-12;
float L=log2(m/0.01); // 0 at |v| = 0.01, one unit per factor of 2
float u=clamp(L/6.0,0.0,1.0);
float iso=L>0.5?1.0-smoothstep(0.7,1.5,abs(fract(L+0.5)-0.5)/max(fwidth(L),1e-6)):0.0; // from |v| = 0.02 up
float fv=fwidth(v);
float node=fv>1e-9?gate*(1.0-smoothstep(0.7,1.5,abs(v)/fv)):0.0;
vec3 c=mix(vec3(0.3,0.32,0.38),base*(0.3+0.9*u),smoothstep(0.0,0.12,u)); // neutral where the field is negligible
c=mix(c,vec3(1.0),0.5*iso);
c=mix(c,vec3(1.0),0.85*node);
float a=max(0.92*(uPlField==0?smoothstep(0.0,0.2,u):smoothstep(0.1,0.5,u)),max(0.6*iso,0.5*node));
return vec4(c,a);
}`;

// The cut drawn flat: uC + x*uE1 + y*uE2 over [-uHalf, uHalf]^2. uRaw = 1 writes the bare
// field value instead of the colour (for checking it against a CPU evaluation).
const sliceSource = (maxAo) => `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
#define MAXAO ${maxAo}
out vec4 o;
uniform vec3 uC,uE1,uE2;
uniform float uHalf,uN;
uniform int uRaw;
${BASIS_GLSL}
${CUT_GLSL}
void main(){
vec2 q=(gl_FragCoord.xy/uN*2.0-1.0)*uHalf;
vec3 p=uC+q.x*uE1+q.y*uE2;
if(uRaw==1){o=vec4(cutField(p),0.0,0.0,1.0);return;}
vec4 c=cutColor(p);
o=vec4(pow(mix(vec3(0.006,0.007,0.012),c.rgb,c.a),vec3(1.0/2.2)),1.0);
}`;

// grid = true: the body density is read from a 3D texture of log(rho) filled once per
// geometry; grid = false: it is evaluated exactly at every ray step (slow, the reference).
// orbGrid = true: the selected orbital and its gradient are read from 3D textures filled once
// per geometry and orbital; orbGrid = false: they are evaluated from the basis set at every
// ray step. That is the costly part of a frame (some thirty exponentials a step for water,
// more with every atom nearby), and what the grid replaces with one or two texture reads.
const fragmentSource = (maxAo, maxAtoms, grid, orbGrid) => `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
precision highp sampler3D;
#define MAXAO ${maxAo}
#define MAXAT ${maxAtoms}
out vec4 o;
uniform vec2 uR;
uniform vec3 uCam,uRt,uUp,uFw; // uCam.z: the camera's distance; the rest: its axes, see cameraBasis
uniform float uRb,uGhost,uBody;
uniform int uSteps;
uniform sampler3D uRho,uPsiA,uPsiB;
// cutting plane: dot(uPl.xyz, p) = uPl.w; uPlMode 0 off, 1 knife (the side the normal points
// to is removed), 2 sheet (fog on both sides); uPlR: radius of the drawn disc
uniform vec4 uPl;
uniform int uPlMode,uFrags; // uFrags: number of fragments in the scene (1 or 2)
uniform float uPlR;
${BASIS_GLSL}
${CUT_GLSL}
float ign(vec2 p){return fract(52.9829189*fract(dot(p,vec2(0.06711056,0.00583715))));}
float erfa(float x){float x2=x*x;return sign(x)*sqrt(1.0-exp(-x2*(1.2732395+0.147*x2)/(1.0+0.147*x2)));}
// mean of a Gaussian band exp(-((l-lvl)/w)^2) between two samples of l, so that a thin
// shell is not missed between ray steps
float bandI(float l0,float l1,float lvl,float w){
float x0=(l0-lvl)/w,x1=(l1-lvl)/w,dx=x1-x0;
if(abs(dx)<1e-3)return exp(-x0*x0);
return 0.8862269*(erfa(x1)-erfa(x0))/dx;
}
float logDensity(vec3 p){
${grid ? 'return texture(uRho,p/uRb*0.5+0.5).x;' : 'return log(density(p)+1e-13);'}
}
void main(){
vec2 uv=(gl_FragCoord.xy*2.0-uR)/uR.y;
vec3 fw=uFw,rt=uRt,up=uUp;
vec3 ro=-uCam.z*fw;
vec3 rd=normalize(fw*2.3+rt*uv.x+up*uv.y);
vec3 L=normalize(rt*0.55+up*0.75-fw*0.45);
vec3 bg=vec3(0.006,0.007,0.012);
float b=dot(ro,rd);
float h=b*b-(dot(ro,ro)-uRb*uRb);
vec3 col=vec3(0.0);
float T=1.0;
float tN[MAXAT];float TN[MAXAT];
for(int k=0;k<uNA;k++){tN[k]=dot(texelFetch(uAtom,ivec2(k,0),0).xyz-ro,rd);TN[k]=1.0;}
float gf=mix(1.0,0.45,uGhost);
// where the ray meets the cutting plane, and what the cut looks like there
float dn=dot(rd,uPl.xyz),s0=dot(ro,uPl.xyz)-uPl.w;
float tp=abs(dn)>1e-6?-s0/dn:-1.0;
vec4 pc=vec4(0.0);
if(uPlMode!=0){
vec3 pp=ro+rd*clamp(tp,0.0,1e4);
pc=cutColor(pp);
float rr=length(pp-uPl.xyz*uPl.w)/uPlR;
float ring=1.0-smoothstep(0.0,0.012,abs(rr-0.985));
pc.rgb=mix(pc.rgb,vec3(0.75,0.8,0.9),ring);
pc.a=max(max(pc.a,0.035),0.3*ring); // the sheet itself is faintly visible, with a rim
if(uPlMode==2)pc.a*=0.62; // "sheet": the fog behind stays visible
if(rr>1.0||tp<=0.0)pc.a=0.0;
}
bool pd=pc.a<=0.0; // cut already composited (or nothing to composite)
if(h>0.0){
h=sqrt(h);
float t0=max(-b-h,0.0),t1=-b+h;
if(uPlMode==1){
// keep only the half-space dot(n,p) <= w
if(dn>1e-6)t1=min(t1,tp);else if(dn<-1e-6)t0=max(t0,tp);else if(s0>0.0)t1=t0;
t1=max(t1,t0);
}
if(!pd&&tp<=t0){col+=T*pc.a*pc.rgb;T*=1.0-pc.a;pd=true;}
float dt=(t1-t0)/float(uSteps);
float t=t0+dt*ign(gl_FragCoord.xy);
float lp=-30.0,lb=-30.0;
for(int i=0;i<uSteps;i++){
vec3 p=ro+rd*t;
float pa,pb;vec3 ga,gb;
${orbGrid ? `vec3 pq=p/uRb*0.5+0.5;
vec4 A=textureLod(uPsiA,pq,0.0);pa=A.x;ga=A.yzw;pb=0.0;gb=vec3(0.0);
if(uFrags>1){vec4 B=textureLod(uPsiB,pq,0.0);pb=B.x;gb=B.yzw;}` : 'orbital2(p,pa,ga,pb,gb);'}
float psi=pa+pb;vec3 g=ga+gb;
float l=log(psi*psi+1e-13),lr=uBody>0.5?logDensity(p):-30.0;
if(i==0){lp=l;lb=lr;}
float so=gf*(3.2*bandI(lp,l,-4.42,0.42)+1.1*bandI(lp,l,-6.44,0.42));
vec3 c=mix(psi>0.0?CPc:CNc,vec3(1.0),0.3*uGhost);
vec3 n=-normalize(psi*g+vec3(1e-12));
float dd=max(dot(n,L),0.0);
float rim=pow(1.0-abs(dot(n,rd)),2.5);
vec3 co=c*(0.28+0.85*dd)+0.3*rim*c;
float sb=uBody>0.5?0.8*bandI(lb,lr,-2.3,0.45):0.0;
vec3 cb=vec3(0.36,0.38,0.43);
// where the parts of the orbital on the two fragments overlap: a third colour, by the
// sign of their product (in phase: bonding, in antiphase: antibonding)
float sv=0.0;vec3 cv=vec3(0.0);
if(uFrags>1){
float ov=pa*pb;
float m=smoothstep(0.0001,0.0015,abs(ov)); // a forming bond gives products of 0.0005–0.02
sv=4.5*m*gf;
vec3 nv=-normalize(sign(ov)*(pa*gb+pb*ga)+vec3(1e-12));
vec3 c3=ov>0.0?vec3(0.35,1.0,0.6):vec3(1.0,0.2,0.65);
cv=c3*(0.28+0.85*max(dot(nv,L),0.0))+0.3*pow(1.0-abs(dot(nv,rd)),2.5)*c3;
so*=1.0-0.8*m;
}
float s=so+sb+sv;
float a=1.0-exp(-s*dt);
col+=T*a*(so*co+sb*cb+sv*cv)/max(s,1e-6);
T*=1.0-a;
lp=l;lb=lr;
t+=dt;
for(int k=0;k<uNA;k++)if(t<tN[k])TN[k]=T;
if(!pd&&t>=tp){col+=T*pc.a*pc.rgb;T*=1.0-pc.a;pd=true;}
if(T<0.01)break;
}
if(!pd){col+=T*pc.a*pc.rgb;T*=1.0-pc.a;for(int k=0;k<uNA;k++)if(tN[k]>tp)TN[k]=T;}
}
col+=T*bg;
float dsz=0.00007*uCam.z*uCam.z;
for(int k=0;k<uNA;k++){
vec3 v=texelFetch(uAtom,ivec2(k,0),0).xyz-ro;
float d2=dot(v,v)-tN[k]*tN[k];
col+=vec3(1.0)*exp(-d2/dsz)*(0.25+0.75*TN[k]);
}
o=vec4(pow(clamp(col,0.0,1.0),vec3(1.0/2.2)),1.0);
}`;

const UNIFORMS = ['uR', 'uCam', 'uRt', 'uUp', 'uFw', 'uRb', 'uGhost', 'uBody', 'uNA', 'uNO', 'uSteps', 'uAtom', 'uShell', 'uPrim', 'uCoef', 'uOcc', 'uRho', 'uPsiA', 'uPsiB', 'uN', 'uCols', 'uPart',
  'uPl', 'uPlMode', 'uFrags', 'uPlField', 'uPlR', 'uC', 'uE1', 'uE2', 'uHalf', 'uRaw'];
const UNITS = { uAtom: 0, uShell: 1, uPrim: 2, uCoef: 3, uOcc: 8 }; // texture units of the data textures
const GRID_SPACING = 0.18; // bohr, target voxel size of the body-density and orbital grids
const bucket = (n, min) => { let b = min; while (b < n) b *= 2; return b; };

// Distance beyond which every function of a shell stays below `tol` in absolute value.
function shellCutoff(shell, tol) {
  const value = (r) => {
    let v = 0;
    for (let k = 0; k < shell.alpha.length; k++) v += Math.abs(shell.coeff[k]) * Math.exp(-shell.alpha[k] * r * r);
    return v * r ** shell.l;
  };
  let r = 30;
  while (r > 0.5 && value(r - 0.25) < tol) r -= 0.25;
  return r;
}

export class FogRenderer {
  // bodyGrid: read the body density from a 3D texture (fast). false = exact evaluation at
  // every ray step, as in the reference demo.
  // orbGrid: read the selected orbital from 3D textures as well, in the frames that are not
  // asked to be exact (see draw). false = every frame is exact.
  constructor(canvas, { bodyGrid = true, orbGrid = true } = {}) {
    this.canvas = canvas;
    const gl = canvas.getContext('webgl2', { antialias: false, alpha: false });
    if (!gl) throw new Error('WebGL2 is not available in this browser.');
    this.gl = gl;
    this.programs = new Map();
    this.textures = {};
    for (const name of Object.keys(UNITS)) {
      const t = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, t);
      for (const [k, v] of [[gl.TEXTURE_MIN_FILTER, gl.NEAREST], [gl.TEXTURE_MAG_FILTER, gl.NEAREST], [gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE], [gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE]])
        gl.texParameteri(gl.TEXTURE_2D, k, v);
      this.textures[name] = t;
    }
    this.scene = null;
    // the grids need rendering into a half-float texture
    const floatTarget = !!(gl.getExtension('EXT_color_buffer_float') || gl.getExtension('EXT_color_buffer_half_float'));
    this.bodyGrid = bodyGrid && floatTarget;
    this.grid = { texture: gl.createTexture(), n: 0 };
    // a: psi and its gradient (the part on fragment 0); b: the part on fragment 1, made only
    // for a scene of two fragments
    this.orbGrid = orbGrid && floatTarget;
    this.orb = { a: null, b: null, n: 0 };
    // the flat sheet the grids are computed on, and the framebuffer their slices are copied through
    this.sheet = { texture: null, fbo: gl.createFramebuffer(), layer: gl.createFramebuffer(), n: 0, cols: 0 };
    this.plane = null;
    this.slice = { fbo: gl.createFramebuffer(), texture: null, key: '' };
    this.canReadFloat = !!gl.getExtension('EXT_color_buffer_float');
  }

  // plane: null, or { normal: [x,y,z] (unit), d (normal . p = d, bohr), mode: 'knife' | 'sheet',
  // field: 'psi' | 'rho' }
  setPlane(plane) {
    this.plane = plane;
  }

  // The cut drawn flat on a size x size image around `center`, spanned by unit vectors e1
  // (right) and e2 (up), half-width `half` bohr. Returns RGBA bytes, top row first; with
  // raw = true a Float32Array of field values instead (bottom row first), or null if the
  // browser cannot read float pixels.
  renderSlice(size, center, e1, e2, half, field, raw = false) {
    const gl = this.gl, s = this.scene, sl = this.slice;
    if (raw && !this.canReadFloat) return null;
    const key = `${size}/${raw}`;
    gl.activeTexture(gl.TEXTURE5);
    if (sl.key !== key) {
      if (sl.texture) gl.deleteTexture(sl.texture);
      sl.texture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, sl.texture);
      gl.texStorage2D(gl.TEXTURE_2D, 1, raw ? gl.R32F : gl.RGBA8, size, size);
      sl.key = key;
    }
    const { prog, U } = this._program(bucket(s.nao, 16), 0, 'slice');
    gl.useProgram(prog);
    this._bind(U);
    gl.uniform3fv(U.uC, center);
    gl.uniform3fv(U.uE1, e1);
    gl.uniform3fv(U.uE2, e2);
    gl.uniform1f(U.uHalf, half);
    gl.uniform1f(U.uN, size);
    gl.uniform1i(U.uPlField, field === 'rho' ? 1 : 0);
    gl.uniform1i(U.uRaw, raw ? 1 : 0);
    gl.bindFramebuffer(gl.FRAMEBUFFER, sl.fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, sl.texture, 0);
    gl.viewport(0, 0, size, size);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    let out;
    if (raw) {
      out = new Float32Array(size * size);
      gl.readPixels(0, 0, size, size, gl.RED, gl.FLOAT, out);
    } else {
      const px = new Uint8Array(4 * size * size);
      gl.readPixels(0, 0, size, size, gl.RGBA, gl.UNSIGNED_BYTE, px);
      out = new Uint8ClampedArray(4 * size * size);
      for (let y = 0; y < size; y++) out.set(px.subarray(4 * size * (size - 1 - y), 4 * size * (size - y)), 4 * size * y);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return out;
  }

  rendererName() {
    const gl = this.gl, ext = gl.getExtension('WEBGL_debug_renderer_info');
    return ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
  }

  // kind: 'grid' (fills the density grid), 'orbgrid' (fills the orbital grid), 'slice' (the
  // flat cut), or the ray march: 'march' followed by 'b' if the body is read from its grid
  // and by 'o' if the orbital is
  _program(maxAo, maxAtoms, kind) {
    const key = `${kind}/${maxAo}/${maxAtoms}`;
    let p = this.programs.get(key);
    if (p) return p;
    const gl = this.gl;
    const shader = (type, src) => {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error('The shader failed to build: ' + gl.getShaderInfoLog(s));
      return s;
    };
    const prog = gl.createProgram();
    gl.attachShader(prog, shader(gl.VERTEX_SHADER, VS));
    gl.attachShader(prog, shader(gl.FRAGMENT_SHADER, kind === 'grid' ? gridSource(maxAo) : kind === 'orbgrid' ? orbitalGridSource(maxAo) : kind === 'slice' ? sliceSource(maxAo)
      : fragmentSource(maxAo, maxAtoms, kind.includes('b', 5), kind.includes('o', 5))));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error('The shader failed to build: ' + gl.getProgramInfoLog(prog));
    p = { prog, U: {} };
    for (const n of UNIFORMS) p.U[n] = gl.getUniformLocation(prog, n);
    this.programs.set(key, p);
    return p;
  }

  // the ray march for the scene: with the orbital from its grid, or exact
  _march(orbGrid) {
    const s = this.scene;
    return this._program(bucket(s.nao, 16), bucket(s.nat, 8), `march${this.bodyGrid ? 'b' : ''}${orbGrid ? 'o' : ''}`);
  }

  _upload(name, width, height, internal, format, data) {
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE0 + UNITS[name]);
    gl.bindTexture(gl.TEXTURE_2D, this.textures[name]);
    gl.texImage2D(gl.TEXTURE_2D, 0, internal, width, height, 0, format, gl.FLOAT, data);
  }

  // positions: bohr, flat 3n, already centred for viewing. result: an engine result.
  // Returns the camera distance and bounding radius chosen for this molecule.
  // During dynamics the caller fixes `extent` (radius that the nuclei never leave) and
  // `planeRadius`, so that the ray step and the density grid do not change from frame to frame.
  // fragments: for every atom 0 or 1, the molecule it belongs to (default: all 0).
  setScene(positions, result, { extent, planeRadius, fragments } = {}) {
    const gl = this.gl, nat = positions.length / 3, { shells, nao, nmo } = result;
    const firstShell = new Int32Array(nat).fill(-1), shellCount = new Int32Array(nat), cut = new Float64Array(nat);
    const shellTex = new Float32Array(4 * shells.length), prims = [];
    let ao = 0;
    shells.forEach((sh, s) => {
      // shells of one atom are contiguous in every engine
      if (firstShell[sh.atom] < 0) firstShell[sh.atom] = s;
      shellCount[sh.atom]++;
      cut[sh.atom] = Math.max(cut[sh.atom], shellCutoff(sh, 3e-4));
      if (sh.l > 2) throw new Error(`shell with l = ${sh.l} is not supported`);
      shellTex.set([sh.l, prims.length / 4, sh.alpha.length, ao], 4 * s);
      for (let k = 0; k < sh.alpha.length; k++) prims.push(sh.alpha[k], sh.coeff[k], 0, 0);
      ao += 2 * sh.l + 1;
    });
    if (ao !== nao) throw new Error(`basis has ${ao} functions, the engine reported ${nao}`);
    const atomTex = new Float32Array(8 * nat);
    let mr = 0;
    for (let k = 0; k < nat; k++) {
      const [x, y, z] = positions.subarray(3 * k, 3 * k + 3);
      atomTex.set([x, y, z, cut[k] * cut[k]], 4 * k);
      atomTex.set([firstShell[k], shellCount[k], fragments ? fragments[k] : 0, 0], 4 * (nat + k));
      mr = Math.max(mr, Math.hypot(x, y, z));
    }
    this._upload('uAtom', nat, 2, gl.RGBA32F, gl.RGBA, atomTex);
    this._upload('uShell', shells.length, 1, gl.RGBA32F, gl.RGBA, shellTex);
    this._upload('uPrim', prims.length / 4, 1, gl.RGBA32F, gl.RGBA, new Float32Array(prims));

    // occupied orbitals, each scaled by sqrt(occupation): rho = sum of squares
    const occupied = [];
    for (let m = 0; m < nmo; m++) if (result.occ[m] > 1e-6) occupied.push(m);
    const occRows = Math.max(1, Math.ceil(occupied.length / 4)), occ = new Float32Array(4 * nao * occRows);
    occupied.forEach((m, k) => {
      const w = Math.sqrt(result.occ[m]), at = 4 * nao * (k >> 2) + (k & 3);
      for (let i = 0; i < nao; i++) occ[at + 4 * i] = w * result.coeff[m * nao + i];
    });
    if (extent) mr = extent;
    const rb = mr + 4.5;
    this.scene = {
      nat, nao, coef: new Float32Array(nao), occRows, result, rb,
      nFragments: fragments && fragments.includes(1) ? 2 : 1,
      // same framing as the reference demo; more ray steps only for larger molecules
      distance: 2.3 * (mr + 2.6),
      planeRadius: planeRadius ?? mr + 2.4, // radius of the drawn cutting disc and half-width of the flat plot
      steps: Math.max(120, Math.min(240, Math.round(120 * rb / 8.5))),
      ghost: 0,
      selected: 0,
      // the grids are filled when a frame first needs them: dynamics may bring several
      // geometries between two frames, and the body may be switched off
      orbStale: true,
      bodyStale: true,
    };
    this._upload('uCoef', nao, 1, gl.R32F, gl.RED, this.scene.coef);
    this._upload('uOcc', nao, occRows, gl.RGBA32F, gl.RGBA, occ);
    return this.scene;
  }

  _bind(U) {
    const gl = this.gl, s = this.scene;
    gl.uniform1f(U.uRb, s.rb);
    gl.uniform1i(U.uNA, s.nat);
    gl.uniform1i(U.uNO, s.occRows);
    Object.entries(UNITS).forEach(([name, unit]) => {
      gl.activeTexture(gl.TEXTURE0 + unit);
      gl.bindTexture(gl.TEXTURE_2D, this.textures[name]);
      gl.uniform1i(U[name], unit);
    });
  }

  // voxels along the edge of the cube around the molecule that the grids cover
  _gridSize() {
    return Math.max(64, Math.min(160, Math.ceil(2 * this.scene.rb / GRID_SPACING)));
  }

  // a cubic 3D texture read with linear interpolation, left bound to the active unit
  _texture3D(format, n) {
    const gl = this.gl, t = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_3D, t);
    gl.texStorage3D(gl.TEXTURE_3D, 1, format, n, n, n);
    for (const [k, v] of [[gl.TEXTURE_MIN_FILTER, gl.LINEAR], [gl.TEXTURE_MAG_FILTER, gl.LINEAR], [gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE], [gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE], [gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE]])
      gl.texParameteri(gl.TEXTURE_3D, k, v);
    return t;
  }

  // Fills the 3D texture `target` (n voxels along the edge) with what the program `kind`
  // computes: one draw call onto the sheet, then a copy of every slice. False if the browser
  // cannot do that after all.
  _bake(kind, target, n, part = 0) {
    const gl = this.gl, s = this.scene, sh = this.sheet;
    if (sh.n !== n) {
      const cols = Math.ceil(Math.sqrt(n)), rows = Math.ceil(n / cols);
      if (cols * n > gl.getParameter(gl.MAX_TEXTURE_SIZE)) return false;
      gl.activeTexture(gl.TEXTURE5);
      gl.deleteTexture(sh.texture);
      sh.texture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, sh.texture);
      gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA16F, cols * n, rows * n);
      gl.bindFramebuffer(gl.FRAMEBUFFER, sh.fbo);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, sh.texture, 0);
      Object.assign(sh, { n, cols, rows });
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, sh.fbo);
    const usable = gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE;
    if (usable) {
      const { prog, U } = this._program(bucket(s.nao, 16), 0, kind);
      gl.useProgram(prog);
      this._bind(U);
      gl.uniform1f(U.uN, n);
      gl.uniform1f(U.uCols, sh.cols);
      gl.uniform1i(U.uPart, part);
      gl.viewport(0, 0, sh.cols * n, sh.rows * n);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, sh.fbo);
      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, sh.layer);
      for (let z = 0; z < n; z++) {
        const x = (z % sh.cols) * n, y = Math.floor(z / sh.cols) * n;
        gl.framebufferTextureLayer(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, target, 0, z);
        gl.blitFramebuffer(x, y, x + n, y + n, 0, 0, n, n, gl.COLOR_BUFFER_BIT, gl.NEAREST);
      }
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return usable && gl.getError() === gl.NO_ERROR;
  }

  // The selected orbital and its gradient on the cube around the molecule. Costs about as
  // much as a hundredth of an exact frame, so it is simply redone for every new geometry and
  // every newly selected orbital.
  _fillOrbital() {
    const gl = this.gl, s = this.scene, o = this.orb, n = this._gridSize(), two = s.nFragments > 1;
    gl.activeTexture(gl.TEXTURE6);
    if (o.n !== n) {
      gl.deleteTexture(o.a);
      gl.deleteTexture(o.b);
      o.a = this._texture3D(gl.RGBA16F, n);
      o.b = null;
      o.n = n;
    }
    if (two && !o.b) o.b = this._texture3D(gl.RGBA16F, n);
    s.orbStale = false;
    // if it cannot be done, every frame is exact from now on
    if (!this._bake('orbgrid', o.a, n, 0) || (two && !this._bake('orbgrid', o.b, n, 1))) this.orbGrid = false;
  }

  // log(rho) on the same cube
  _fillGrid() {
    const gl = this.gl, s = this.scene, g = this.grid, n = this._gridSize();
    gl.activeTexture(gl.TEXTURE4);
    if (g.n !== n) {
      gl.deleteTexture(g.texture);
      g.texture = this._texture3D(gl.R16F, n);
      g.n = n;
    }
    s.bodyStale = false;
    // if it cannot be done, the body is evaluated exactly from now on
    if (!this._bake('grid', g.texture, n)) this.bodyGrid = false;
  }

  // Show orbital m in colour.
  select(m) {
    const s = this.scene, gl = this.gl, { result } = s;
    s.selected = m;
    s.coef.set(result.coeff.subarray(m * s.nao, (m + 1) * s.nao), 0);
    s.ghost = 1 - Math.min(1, result.occ[m] / 2); // empty orbitals are drawn paler
    this._upload('uCoef', s.nao, 1, gl.R32F, gl.RED, s.coef);
    s.orbStale = true;
  }

  // the empty background, for a scene that has no cloud to show
  clear() {
    const gl = this.gl;
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    gl.clearColor(0.098, 0.105, 0.134, 1); // the shader's background after its gamma
    gl.clear(gl.COLOR_BUFFER_BIT);
  }

  // camera: { yaw, pitch, distance } or one turned freely, see turnCamera
  // exact: evaluate the orbital from the basis set at every ray step instead of reading its
  // grid. Many times slower; meant for the frame that stays on the screen.
  draw(camera, { body = true, exact = false } = {}) {
    const s = this.scene, gl = this.gl;
    if (!s) return;
    // either may find its grid unusable
    if (this.bodyGrid && body && s.bodyStale) this._fillGrid();
    if (this.orbGrid && !exact && s.orbStale) this._fillOrbital();
    const fromGrid = this.orbGrid && !exact;
    const { prog, U } = this._march(fromGrid);
    gl.useProgram(prog);
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    gl.uniform2f(U.uR, this.canvas.width, this.canvas.height);
    const basis = cameraBasis(camera);
    gl.uniform3f(U.uCam, 0, 0, camera.distance);
    gl.uniform3fv(U.uRt, basis.rt);
    gl.uniform3fv(U.uUp, basis.up);
    gl.uniform3fv(U.uFw, basis.fw);
    this._bind(U);
    gl.uniform1f(U.uGhost, s.ghost);
    gl.uniform1f(U.uBody, body ? 1 : 0);
    gl.uniform1i(U.uSteps, s.steps);
    gl.activeTexture(gl.TEXTURE4);
    gl.bindTexture(gl.TEXTURE_3D, this.grid.texture);
    gl.uniform1i(U.uRho, 4);
    if (fromGrid) {
      gl.activeTexture(gl.TEXTURE6);
      gl.bindTexture(gl.TEXTURE_3D, this.orb.a);
      gl.uniform1i(U.uPsiA, 6);
      gl.activeTexture(gl.TEXTURE7);
      gl.bindTexture(gl.TEXTURE_3D, s.nFragments > 1 ? this.orb.b : this.orb.a);
      gl.uniform1i(U.uPsiB, 7);
    }
    const pl = this.plane;
    gl.uniform1i(U.uFrags, s.nFragments);
    gl.uniform1i(U.uPlMode, !pl ? 0 : pl.mode === 'knife' ? 1 : 2);
    if (pl) {
      gl.uniform4f(U.uPl, pl.normal[0], pl.normal[1], pl.normal[2], pl.d);
      gl.uniform1i(U.uPlField, pl.field === 'rho' ? 1 : 0);
      gl.uniform1f(U.uPlR, s.planeRadius);
    }
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }
}

// Camera position and axes; the shader receives the axes from here. A camera is
// { yaw, pitch, distance }, looking at the origin with the y axis up, until it has been turned
// by hand: then it carries its own axes (see turnCamera) and they are what counts.
export function cameraBasis(camera) {
  if (camera.axes) { const { fw, rt, up } = camera.axes; return { ro: fw.map((v) => -v * camera.distance), fw, rt, up }; }
  const { yaw, pitch, distance } = camera, cp = Math.cos(pitch);
  const ro = [distance * Math.sin(yaw) * cp, distance * Math.sin(pitch), distance * Math.cos(yaw) * cp];
  const fw = ro.map((v) => -v / distance);
  const rl = Math.hypot(fw[2], fw[0]), rt = [-fw[2] / rl, 0, fw[0] / rl];
  const up = [rt[1] * fw[2] - rt[2] * fw[1], rt[2] * fw[0] - rt[0] * fw[2], rt[0] * fw[1] - rt[1] * fw[0]];
  return { ro, fw, rt, up };
}

// The camera turned about the origin: by `aboutUp` radians about its own up axis, then by
// `aboutRight` about its own right axis. Being its own axes, they have no pole to stop at: the
// view can be turned on and on in any direction. yaw and pitch are kept in step for reference.
export function turnCamera(camera, aboutUp, aboutRight) {
  let { fw, rt, up } = cameraBasis(camera);
  const turn = (v, axis, angle) => {
    const cs = Math.cos(angle), sn = Math.sin(angle), d = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2];
    const cr = [axis[1] * v[2] - axis[2] * v[1], axis[2] * v[0] - axis[0] * v[2], axis[0] * v[1] - axis[1] * v[0]];
    return v.map((x, i) => x * cs + cr[i] * sn + axis[i] * d * (1 - cs));
  };
  fw = turn(fw, up, aboutUp); rt = turn(rt, up, aboutUp);
  fw = turn(fw, rt, aboutRight);
  // keep the axes exactly orthonormal over thousands of small turns
  const fl = Math.hypot(...fw);
  fw = fw.map((v) => v / fl);
  const along = rt[0] * fw[0] + rt[1] * fw[1] + rt[2] * fw[2];
  rt = rt.map((v, i) => v - along * fw[i]);
  const rl = Math.hypot(...rt);
  rt = rt.map((v) => v / rl);
  up = [rt[1] * fw[2] - rt[2] * fw[1], rt[2] * fw[0] - rt[0] * fw[2], rt[0] * fw[1] - rt[1] * fw[0]];
  return { ...camera, axes: { fw, rt, up }, yaw: Math.atan2(-fw[0], -fw[2]), pitch: Math.asin(Math.max(-1, Math.min(1, -fw[1]))) };
}

// Screen position (CSS pixels) of a point for the same camera as the shader.
export function project(point, camera, width, height) {
  const { ro, fw, rt, up } = cameraBasis(camera);
  const v = [point[0] - ro[0], point[1] - ro[1], point[2] - ro[2]];
  const z = v[0] * fw[0] + v[1] * fw[1] + v[2] * fw[2];
  const x = (v[0] * rt[0] + v[1] * rt[1] + v[2] * rt[2]) / z * 2.3, y = (v[0] * up[0] + v[1] * up[1] + v[2] * up[2]) / z * 2.3;
  return [width / 2 + x * height / 2, height / 2 - y * height / 2, z];
}

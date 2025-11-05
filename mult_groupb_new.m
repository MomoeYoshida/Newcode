%
% do a multi-plane / multi-measurement example
% here, specifying all overlap parameters in detail by user, not inferred
%
% Scale parameter:
%   [sc, rowinfo, colinfo]
%   info:  [bandsize numbands interoverlap bdyoverlap startindex wrap]
%
%   Rules:  startindex <= 1
%           if wrap is on (1) - bdyoverlap is even, startind = 1-bdyov/2
%                               numbands*(bandsiz-interov)+interov must be even
%

function [r1,r2,r3,r4,r5,r6] = mult_groupb_new( sc, st_dense, varargin );

if (length(sc)==1),
  % Regular call to mult_fields6, change nothing
  comstr = '';
  if (nargout > 0),
    comstr='[r1';
    for i=2:nargout, comstr=[comstr ',r' int2str(i)]; end;
    comstr=[comstr '] '];
  end;
%  comstr = [comstr '= mult_fields6_new( sc, st_dense'];
  comstr = [comstr '= mult_fields6_newest( sc, st_dense'];
  for i=1:(nargin-2), comstr=[comstr ',varargin{' int2str(i) '}']; end;
  comstr = [comstr ');']
  
  % Do the call
  eval(comstr);
  return;
end;

% We're requesting a decomposition, will need to break down matrices

% Get parameter positions
i_nf = 1;
num_fields = eval(['varargin{' num2str(i_nf) '}(1)']);
i_nw = i_nf+1+4*num_fields;
i_mm = i_nw+1+eval(['varargin{' num2str(i_nw) '}']);
num_meas = size(eval(['varargin{' num2str(i_mm) '}']),1);
i_re = i_mm+1+2*num_meas;
eval(['opts = varargin{' num2str(i_re+2) '};']);
meas_siz = eval(['size(varargin{' num2str(i_mm+1) '})']);

% Decode decomposition request
if length(sc)<13
  error( 'Too few parameters in sc.' );
end;

sizr = sc(2); sizc = sc(8);
numr = sc(3); numc = sc(9);
or   = sc(4); oc   = sc(10);
bor  = sc(5); boc  = sc(11);
ofsr = sc(6); ofsc = sc(12);
wrr  = sc(7); wrc  = sc(13);

if (ofsr > 1 | ofsc > 1),
  error( 'Must have offsets at left edge or earlier.' );
end;
if (wrr > 0.5),
  if (rem(bor,2)~=0), error('Row boundary overlap must be even.'); end;
  if (ofsr ~= 1-bor/2), error('Wrong row initial offset.'); end;
  if (rem(numr*(sizr-or)+or,2)~=0), error('Row overlapped domain size is odd.' ); end;
end;
if (wrc > 0.5),
  if (rem(boc,2)~=0), error('Column boundary overlap must be even.'); end;
  if (ofsc ~= 1-boc/2), error('Wrong column initial offset.'); end;
  if (rem(numc*(sizc-oc)+oc,2)~=0), error('Column overlapped domain size is odd.' ); end;
end;
exr = numr*(sizr-or)+or-meas_siz(1);
if (wrr < 0.5),
  if (exr < 0), error( 'Row bands do not cover domain.' ); end;
else,
  if (exr ~= bor), error( 'Wrapped row bands do not fit domain.' ); end;
end;

exc = numc*(sizc-oc)+oc-meas_siz(2);
if (wrc < 0.5),
  if (exc < 0), error( 'Column bands do not cover domain.' ); end;
else,
  if (exc ~= boc), error( 'Wrapped column bands do not fit domain.' ); end;
end;


% Allocate collage matrices
for i=1:nargout,
  eval(['cout' int2str(i) ' = zeros(sizr*numr,sizc*numc);']);
end;

% Loop over problem
for r=1:numr, for c=1:numc,
disp([r c])
  comstr = '[ r1';
  for p=2:nargout, comstr=[comstr ',r' int2str(p)]; end;
%  comstr = [comstr '] = mult_fields6_new( sc(1), st_dense'];
  comstr = [comstr '] = mult_fields6_newest( sc(1), st_dense'];

  % translate variables
  for p=1:(nargin-2),
    if rem(eval(['size(varargin{' int2str(p) '})']),meas_siz),
      comstr = [comstr ',varargin{' int2str(p) '}'];
    else,
      comstr = [comstr ',sub' int2str(p)];
      rr = ofsr+(0:(sizr-1))+(sizr-or)*(r-1);
      rr(find(rr<1)) = rr(find(rr<1)) + meas_siz(1);
      rr(find(rr>meas_siz(1))) = rr(find(rr>meas_siz(1))) - meas_siz(1);
      rc = ofsc+(0:(sizc-1))+(sizc-oc)*(c-1);
      rc(find(rc<1)) = rc(find(rc<1)) + meas_siz(2);
      rc(find(rc>meas_siz(2))) = rc(find(rc>meas_siz(2))) - meas_siz(2);
      if (eval(['size(varargin{' int2str(p) '},2)']) == 3*meas_siz(2)),
        eval(['sub' int2str(p) ' = [varargin{' int2str(p) '}(rr,rc) varargin{' int2str(p) '}(rr,meas_siz(2)+rc) varargin{' int2str(p) '}(rr,meas_siz(2)*2+rc)];']);
      else, 
        eval(['sub' int2str(p) ' = varargin{' int2str(p) '}(rr,rc);']);
      end;
    end;
  end;

  % call routine
  % should really modify options vector to not allocate mem after first run
  comstr = [comstr ');'];
  comstr

  eval(comstr);

%if (min(min(r4))<0),
%disp('*****************************');
%x25=[1 -1];
%eval(comstr);
%x25=1;
%keyboard
%end;

%r1=sub16;r2=sub17;r3=sub18;r4=sub19;r5=sub20;r6=sub21;
%minmax(rr)
%minmax(rc)
%imagesc(r6); pause(0);
%if (min(min([r2 r4 r6]))<0),
%keyboard
%end;

  % copy results into collage
  for p=1:nargout,
    eval( ['cout' int2str(p) '((1:sizr)+(r-1)*sizr,(1:sizc)+(c-1)*sizc) = r' int2str(p) ';'] );
  end;
end; 
end;

save allstuff

% Finally, unoverlap collage
for i=1:nargout,
%  eval(['r' int2str(i) ' = group_overlap( numr, numc, or, oc, cout' int2str(i) ', wrr * bor, wrc * boc );']);
  eval(['r' int2str(i) ' = group_overlap( numr, numc, or, oc, cout' int2str(i) ' );']);
end;


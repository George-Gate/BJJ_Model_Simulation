%Author: George-Gate
%Date: 2015/10/6
%--------------------------------------------------------------------------
% Display a fock state in the form |n1,n2>.
%
% str = displayFockState( k2nn, psi ,tag)
%
%[Input]
%  k2nn - Just pass this array.
%  psi  - The state vector you want to display
%  tag  - A string that will show before the vector.
%[Return]
%  str - A string which is printed to Command Window
%
function str = displayFockState( k2nn, psi ,tag)

    if (nargin==2)
        tag='';
    end
    flag=true;
    str='';
    len=length(psi);
    tol=1e-4*max(abs(psi));
    for i=1:len
        if ( abs(psi(i))>tol )
            if (flag)
                flag=false;
            else
                str=[str,' + '];
            end
            if (abs(imag(psi(i))) > 1e-3/sqrt(len))
                str=[str,'(',num2str(psi(i),4),')'];
            else
                str=[str,num2str(real(psi(i)),4)];
            end
            str=[str,'*|',num2str(k2nn(1,i)),',',num2str(k2nn(2,i)),'>'];
        end
    end
    if flag
        str='0';
    end
    display([tag,str]);
end


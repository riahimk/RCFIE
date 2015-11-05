function []=CREATE_TABULAR(ScPb,tab)
file_name=sprintf('tabular%s.tex',ScPb.type(1));
FID = fopen(file_name, 'a');
if(ScPb.k==1 & ScPb.n==8)
    fprintf(FID, '\\begin{tabular}{|cc|cc|cc|}\\hline \n');
    fprintf(FID, ' & & $ICFIE-\\mathcal{I}$& &$ICFIE-\\mathcal{R}$ & \\\\');
    fprintf(FID, 'k &  n & $\\Re{e}\\left(u_\\infty(d)\\right)$ & $\\Im{m}\\left((u_\\infty)(d)\\right)$ & $\\Re{e}\\left(u_\\infty(d)\\right)$ & $\\Im{m}\\left((u_\\infty)(d)\\right)$ \\\\ \\hline \n');
end

    fprintf(FID, '%d & %d & %f & %f & %f & %f \\\\ \n', ScPb.k, ScPb.n,real(tab(1)),imag(tab(1)),real(tab(2)),imag(tab(2)));
    fprintf(FID, '\\hline ');
    fprintf(FID, '\n');
    if(ScPb.k == 5 &  ScPb.n==128)
    fprintf(FID, '\\end{tabular}\n');
    end
    fclose(FID);
end
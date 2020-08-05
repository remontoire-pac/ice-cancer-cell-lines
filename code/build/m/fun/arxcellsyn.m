function [cn,qd] = arxcellsyn(incn)
% function [cn,qd] = arxcellsyn(incn)
%  Convert cell-line synonyms incn to ArxSpan names cn with max qd entries.
%   Created: PAC 2019-05-06
%   Modified: PAC/SM 2020-03-10
%   Modified: PAC 2020-04-25 to edit documentation.
%   Last updated: PAC 2020-07-06 to deal with multiple 'ArxSpan id' per CElL_SAMPLE_ID.
%

    s = 'cell-line-synon-names.csv';
    T = chkfile(s,'CDDB synonym',',','string',true);
    % find all available ArxSpan id
    q1u = unique(T.CELL_SAMPLE_NAME(contains(T.CELL_NAME_TYPE,'Arx')));
    qd = max(double(strrep(unique(q1u),'ACH-',''))); % max CCL dimension
    % insist input incn is a vector of strings
    assert(isstring(incn),'Input list incn must be have string data type.');
    assert(isvector(incn),'Input list incn must be a vector of strings.');
    ncn = numel(incn);    
    for qi=1:ncn
        qq = unique(T.CELL_SAMPLE_ID(ismember(T.CELL_SAMPLE_NAME,incn(qi))));
        if (numel(qq)>1)
            incn(qi) = "Multiple CELL_SAMPLE_ID found for this CELL_SAMPLE_NAME.";
        elseif (numel(qq)==1)
            qa = unique(T.CELL_SAMPLE_NAME(T.CELL_SAMPLE_ID==qq&T.CELL_NAME_TYPE=="ArxSpan id"));
            if (numel(qa)>1)
                % trap correct ID; brittle: need 'ArxSpan id (pending)' in CDDB.
                fprintf('More than one ArxSpan id for CELL_SAMPLE_ID %i, from name %s.\n',qq,incn(qi));
                incn(qi) = qa(1);
                fprintf('Keeping lower numeric ACH value %s (OK as of 20Q2).\n',incn(qi));
            elseif (numel(qa)==1)
                incn(qi) = qa;
            else
                incn(qi) = "No ArxSpan identifier found for this CELL_SAMPLE_ID.";
            end
        else
            incn(qi) = "No CELL_SAMPLE_ID found for this CELL_SAMPLE_NAME.";
        end
    end
    cn = incn;
    
end

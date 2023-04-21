function ExtremeFinalMatrix = LinkIndependentFormulation(new,Total_Links,Length,Angles)

FinalMatrix = zeros(2,new);
ExtremeFinalMatrix = zeros(2,Total_Links);
syms AdditionMatrix [2 Total_Links]
AdditionMatrix(1:2,1:end) = 0;

    for link_num = 1 : new
    
        CurrentMatrix = Runner(link_num,new,Length,Angles);
        FinalMatrix = FinalMatrix + CurrentMatrix;
        [~,col] = size(FinalMatrix);
    
    end

AdditionMatrix(1,1:col) = FinalMatrix(1,1:end);
AdditionMatrix(2,1:col) = FinalMatrix(2,1:end);

ExtremeFinalMatrix = ExtremeFinalMatrix + AdditionMatrix;

end


function Body = Runner(link_num,Total_Links,Length,Angles)

EndMatrix = zeros(2,Total_Links);

    if link_num ~= Total_Links
        Matrix = JacobianLinks(link_num,Total_Links,Length,Angles);
        Body = EndMatrix + Matrix;
    else
        Matrix = JacobianLinkEnd(link_num,Total_Links,Length,Angles);
        Body = EndMatrix + Matrix;
    end

end 

%%
function [Jac_Link] = JacobianLinkEnd(link_num,Total_Links,Length,Angles)

% syms Angles Length [1 Total_Links]

J_S = zeros(2,Total_Links); % Complete Jacobian Matrix with External Forces
J_L = sym('J_L',[2,Total_Links]); % Jacobian for current link

Total_Angle = sum(Angles);
J_L(1,1:end) = -0.5*Length(link_num)*sin(Total_Angle);
J_L(2,1:end) = 0.5*Length(link_num)*cos(Total_Angle);
Jac_Link = J_S + J_L;

end

%%
function [Jac_Link] = JacobianLinks(link_num,Total_Links,Length,Angles)

% syms Angles Length [1 Total_Links]

J_S = zeros(2,Total_Links); % Complete Jacobian Matrix with External Forces
J_L = sym('J_L',[2,Total_Links]); % Jacobian for current link

Total_Angle = sum(Angles(1:link_num));
J_L(1,1:link_num) = -Length(link_num)*sin(Total_Angle);
J_L(1,link_num+1:end) = 0;
J_L(2,1:link_num) = Length(link_num)*cos(Total_Angle);
J_L(2,link_num+1:end) = 0;
Jac_Link = J_S + J_L;

end

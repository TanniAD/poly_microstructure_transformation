function[A,D] = BoundaryReduction(B)
faces = [1 length(B); 1 length(B);1 length(B)];
C = cell(6, 6);

    C{1,1} = slice(B,[],[],1); C{1,2} = get([C{1,1}],'CData'); % Top
    C{2,1} = slice(B,[],[],length(B)); C{2,2} = get([C{2,1}],'CData'); % Bottom
        C{3,1} = slice(B,1,[],[]); C{3,2} = get([C{3,1}],'CData'); % Left
        C{4,1} = slice(B,length(B),[],[]); C{4,2} = get([C{4,1}],'CData');% Right
    C{5,1} = slice(B,[],1,[]); C{5,2} = get([C{5,1}],'CData'); % Front
    C{6,1} = slice(B,[],length(B),[]); C{6,2} = get([C{6,1}],'CData'); % Back
    
for i = 1:numel(faces)
    [C{i,3},C{i,4}] = count_unique(C{i,2});
    for j = 1:length(C{i,3})
        B(B==C{i,3}(j)) = NaN;
        C{i,5} = B;
    end
end
        A = B;
        [P,Q] = count_unique(A);
        D = nthroot(6*Q/pi,3);
% for k = 1:700
%     figure(1),
%     s(k) = surf(squeeze(A(k,:,:)));
%     view(0,90),shading interp
%     pause(0.1)
% end
end
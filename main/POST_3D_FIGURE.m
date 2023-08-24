function POST_3D_FIGURE(X1, X2, Y1, Y2, step, U, NODE)
    %%%后处理函数，输入xy方向的计算节点，求解出的函数值，来做图像拟合
    NX=(X2-X1)/step;
    NY=(Y2-Y1)/step;
    [X,Y]=meshgrid(X1:NX:X2,Y1:NY:Y2);
    for i=1:size(X,1)
        for j=1:size(Y,2)
            x=X(i,j);
            y=Y(i,j);
            for k=1:size(NODE,1)
                if x==NODE(k,1)&&y==NODE(k,2)
                    Z(i,j)=U(k);
                end
            end
        end
    end
    surf(X,Y,Z);
    grid on 
    shading interp
end
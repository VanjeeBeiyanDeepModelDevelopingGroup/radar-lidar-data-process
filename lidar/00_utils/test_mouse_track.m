% 主函数

function test_mouse_track()
figure;
axis([-10,10,0,5]);
set(gcf,'WindowButtonDownFcn',@ButttonDownFcn);

 

% 回调函数

function ButttonDownFcn(src,event)
pt = get(gca,'CurrentPoint');
x = pt(1,1);
y = pt(1,2);
fprintf('x=%f,y=%f\n',x,y);
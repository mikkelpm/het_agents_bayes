ll93 = load('loglike_space_930.mat');
ll96 = load('loglike_space_960.mat');

F1 = figure;
scatter(ll93.all_loglikes(:),ll96.all_loglikes(:),'.')
refline([1,0])
xlabel('\beta=.93')
ylabel('\beta=.96')
print(F1,'ll93_96','-dpng');

scatter(reshape(repmat(ll93.data_hh(:,:,2),[1,1,500]),[],1),ll93.all_loglikes(:),'.')
xlabel('net income')
ylabel('\beta=.93')
print(F1,'ll93','-dpng');

scatter(reshape(repmat(ll93.data_hh(:,:,2),[1,1,500]),[],1),ll96.all_loglikes(:),'.')
xlabel('net income')
ylabel('\beta=.96')
print(F1,'ll96','-dpng');

scatter(reshape(repmat(ll93.data_hh(:,:,3),[1,1,500]),[],1),ll93.all_loglikes(:),'.')
xlabel('asset')
ylabel('\beta=.93')
print(F1,'ll93_asset','-dpng');

scatter(reshape(repmat(ll93.data_hh(:,:,3),[1,1,500]),[],1),ll96.all_loglikes(:),'.')
xlabel('asset')
ylabel('\beta=.96')
print(F1,'ll96_asset','-dpng');
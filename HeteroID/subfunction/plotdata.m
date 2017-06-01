function plotdata(result, y, A, parameter,figurename)


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%       NON GROUP           %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(parameter.group_onoff, 'off')
    train_length = round(parameter.prop*parameter.observation_dim);
    
    plot_onoff = parameter.plot_onoff;
    diff_onoff = parameter.diff_onoff;
    group_type = parameter.group_type;
    
    if strcmp(group_type,'MMV')
        [p,L] = size(y);
        y = vec(y');
        A = kron(A,eye(L));
    end
    
    
    
    [m,n] = size(A);
    t = [1:m]';
    if strcmp(parameter.group_onoff,'on')
        partition = result.partition;
    end
    true = result.Compare(:,1);
    est_1 = result.Compare(:,2);
    est_end = result.Compare(:,end);
    
    
    if strcmp(diff_onoff,'off')
        
        f = figure('visible',plot_onoff);
        subplot(2,2,1)
        plot(t(1:train_length),y(1:train_length),'r *-','LineWidth',1);hold on;
        plot(t(train_length+1:end),y(train_length+1:end),'g *-','LineWidth',1);hold on;
        plot(t,A*est_1,'b o-','LineWidth', 1.5);hold on;
        legend('oringnal train set','oringnal test set', 'estimated data');
        title('Data: First iteration');xlabel('Time');ylabel('Output');
        line([t(train_length+1) t(train_length+1)], [1.01*min(min(y)) 1.01*max(max(y))],...
            'LineWidth',2,'Color',[.8 .8 .8]);
        
        subplot(2,2,3)
        plot(t(1:train_length),y(1:train_length),'r *-','LineWidth',1);hold on;
        plot(t(train_length+1:end),y(train_length+1:end),'g *-','LineWidth',1);hold on;
        plot(t,A*est_end,'b o-','LineWidth', 1.5);hold on;
        legend('oringnal train set','oringnal test set', 'estimated data');
        title('Data: End iteration');xlabel('Time');ylabel('Output');
        line([t(train_length+1) t(train_length+1)], [1.01*min(min(y)) 1.01*max(max(y))],...
            'LineWidth',2,'Color',[.8 .8 .8]);
        
        subplot(2,2,2)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_1,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('First iteration')
        subplot(2,2,4)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_end,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('End iteration')
        
%         saveas(f, figurename);
        
    elseif strcmp(diff_onoff,'on')
        %%
        y_d = diff(y);
        A_d = diff(A);
        
        %  Plotting
        %     [~,threshold_true] = prunedic(A, y, result.x_true,0);
        %     [~,threshold_estimate] = prunedic(A, y, result.estimate,0);
        %     result.threshold = [threshold_true, threshold_estimate];
        %     result.nmse_y = norm(y - A*result.estimate,2)^2/norm(y,2)^2;
        %     result = orderfields(result);
        
        % save('RESULT.mat');
        % load('RESULT.mat');
        
        if length(y)>=15 && adftest(y) == 1
            f = figure('visible',plot_onoff);
            subplot(2,2,1)
            plot(t,y,'m o-','LineWidth',1.5);hold on;
            plot(t,A*est_1,'b *-','LineWidth', 1.5);hold on;
            legend('oringnal','estimated oringnal');
            title('Result');xlabel('Time');ylabel('Output');
            subplot(2,2,3)
            plot(t,y,'m o-','LineWidth',1.5);hold on;
            plot(t,A*est_end,'b *-','LineWidth', 1.5);hold on;
            legend('oringnal','estimated oringnal');
            title('Result');xlabel('Time');ylabel('Output');
        else
            f = figure('visible',plot_onoff);
            subplot(2,2,1)
            plot(t,y,'m o-','LineWidth',1);hold on;
            plot(t(1:end-1),y_d,'r o-','LineWidth',1.5); hold on;
            %     plot(t,cumsum([y(1); A_d*est_end]),'b-','LineWidth', 2);hold on;
            plot(t,A*est_1,'b *-','LineWidth', 1.5);hold on;
            plot(t(1:end-1), A_d*est_1,'g-','LineWidth', 2);hold on;
            legend('oringnal','diff','estimated oringnal','estimated diff');
            title('Result');xlabel('Time');ylabel('Output');
            subplot(2,2,3)
            plot(t,y,'m o-','LineWidth',1);hold on;
            plot(t(1:end-1),y_d,'r o-','LineWidth',1.5); hold on;
            %     plot(t,cumsum([y(1); A_d*est_end]),'b-','LineWidth', 2);hold on;
            plot(t,A*est_end,'b *-','LineWidth', 1.5);hold on;
            plot(t(1:end-1), A_d*est_end,'g-','LineWidth', 2);hold on;
            legend('oringnal','diff','estimated oringnal','estimated diff');
            title('Result');xlabel('Time');ylabel('Output');
        end
        % display(result);
        
        
        subplot(2,2,2)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_1,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('First iteration')
        
        
        subplot(2,2,4)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_end,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('End iteration')
        
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                 GROUP           %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(parameter.group_onoff, 'on')
    plot_onoff = parameter.plot_onoff;
    diff_onoff = parameter.diff_onoff;
    group_type = parameter.group_type;
    
    if strcmp(group_type,'MMV')
        [p,L] = size(y);
        y = vec(y');
        A = kron(A,eye(L));
    end
    
    
    
    [m,n] = size(A);
    t = [1:m]';
    if strcmp(parameter.group_onoff,'on')
        partition = result.partition;
    end
    true = result.Compare(:,1);
    est_1 = result.Compare(:,2);
    est_end = result.Compare(:,end);
    
    if strcmp(diff_onoff,'on')
        
        y_d = diff(y);
        A_d = diff(A);
         
        %% Plotting
        %     [~,threshold_true] = prunedic(A, y, result.x_true,0);
        %     [~,threshold_estimate] = prunedic(A, y, result.estimate,0);
        %     result.threshold = [threshold_true, threshold_estimate];
        %     result.nmse_y = norm(y - A*result.estimate,2)^2/norm(y,2)^2;
        %     result = orderfields(result);
        
        % save('RESULT.mat');
        % load('RESULT.mat');
        
        if length(y)>=15 && adftest(y) == 1
            f = figure('visible',plot_onoff);
            subplot(2,2,1)
            plot(t,y,'m o-','LineWidth',1.5);hold on;
            plot(t,A*est_1,'b *-','LineWidth', 1.5);hold on;
            legend('oringnal','estimated oringnal');
            title('Result');xlabel('Time');ylabel('Output');
            subplot(2,2,3)
            plot(t,y,'m o-','LineWidth',1.5);hold on;
            plot(t,A*est_end,'b *-','LineWidth', 1.5);hold on;
            legend('oringnal','estimated oringnal');
            title('Result');xlabel('Time');ylabel('Output');
        else
            f = figure('visible',plot_onoff);
            subplot(2,2,1)
            plot(t,y,'m o-','LineWidth',1);hold on;
            plot(t(1:end-1),y_d,'r o-','LineWidth',1.5); hold on;
            %     plot(t,cumsum([y(1); A_d*est_end]),'b-','LineWidth', 2);hold on;
            plot(t,A*est_1,'b *-','LineWidth', 1.5);hold on;
            plot(t(1:end-1), A_d*est_1,'g-','LineWidth', 2);hold on;
            legend('oringnal','diff','estimated oringnal','estimated diff');
            title('Result');xlabel('Time');ylabel('Output');
            subplot(2,2,3)
            plot(t,y,'m o-','LineWidth',1);hold on;
            plot(t(1:end-1),y_d,'r o-','LineWidth',1.5); hold on;
            %     plot(t,cumsum([y(1); A_d*est_end]),'b-','LineWidth', 2);hold on;
            plot(t,A*est_end,'b *-','LineWidth', 1.5);hold on;
            plot(t(1:end-1), A_d*est_end,'g-','LineWidth', 2);hold on;
            legend('oringnal','diff','estimated oringnal','estimated diff');
            title('Result');xlabel('Time');ylabel('Output');
        end
        % display(result);
        
        
        subplot(2,2,2)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_1,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('First iteration')
        
        
        subplot(2,2,4)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_end,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('End iteration')
        
        
    elseif strcmp(diff_onoff,'off')
        
        f = figure('visible',plot_onoff);
        subplot(2,2,1)
        plot(t,y,'m o-','LineWidth',1.5);hold on;
        plot(t,A*est_1,'b o-','LineWidth', 1.5);hold on;
        legend('oringnal','estimated oringnal');
        title('Data: First iteration (Group Lasso)');xlabel('Time');ylabel('Output');
        subplot(2,2,3)
        plot(t,y,'m o-','LineWidth',1.5);hold on;
        plot(t,A*est_end,'b o-','LineWidth', 1.5);hold on;
        legend('oringnal','estimated oringnal');
        title('Data: End iteration (Our Algorithm)');xlabel('Time');ylabel('Output');
        
        subplot(2,2,2)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_1,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('Identification: First iteration (Group Lasso)')
        subplot(2,2,4)
        plot( true,'r','LineWidth', 1.5);
        hold on;
        plot(est_end,'b o-','LineWidth', 1.5)
        hold on;
        legend('true','estimate');
        if strcmp(parameter.group_onoff,'on')
            cumsumpartition = cumsum(partition);
            if strcmp(parameter.group_onoff, 'on')
                for i = 1:length(partition)
                    line([cumsumpartition(i)+0.5, cumsumpartition(i)+0.5], ...
                        [min(min(est_1, true))-0.5, max(max(est_1, true))+0.5]); hold on;
                end
            end
            axis([0.5 length(true)+0.5 min(min(est_1, true))-0.5, max(max(est_1, true))+0.5])
        end
        title('Identification: End iteration (Our Algorithm)')
        
%         saveas(f, figurename);
    end
    
    
    
end



function y = overlap_keep_convolution(x_e, x_c)
    % x_e: 输入信号 (512 点随机序列)
    % x_c: 卷积核 (长度为 8 的随机序列)

    % 获取输入信号和卷积核的长度
    N = length(x_e);
    M = length(x_c);
    
    % 将输入信号分成 8 段
    num_segments = 8;
    segment_length = N / num_segments;
    
    % 补零卷积核到每一段卷积的长度
    x_c_padded = [x_c; zeros(segment_length-1,1)];
    
    % 初始化输出信号 y
    y = zeros(N + M - 1, 1); % 卷积结果的长度

    for i = 1:num_segments+1
        % 提取当前段并补零到卷积核长度
        x_segment = x_e((i-1)*segment_length+1 : min(i*segment_length, N));
        if i==1
            x_segment_padded =[zeros(M-1,1);x_segment];%前面补7个0
        end
        if i>1
            x_keep_segment=x_e((i-1)*segment_length-6:(i-1)*segment_length);
            x_segment_padded =[x_keep_segment;x_segment];%前面保留7个数
        end
         if i==num_segments+1%增加一段全零样本
            x_keep_segment=x_e((i-1)*segment_length-6:(i-1)*segment_length);
            x_segment_padded =[x_keep_segment;zeros(M+N/8-8,1)];
        end
        
        % 使用FFT计算卷积
        X_segment = (fft(x_segment_padded));
        X_c = (fft(x_c_padded));
        convolution_result = (ifft(X_segment.*X_c));
        
        % 把每一段的前七个删除然后拼接
        if i<=num_segments
        output_start =(i-1)*segment_length+1;
        output_end = output_start + segment_length-1;
        y(output_start : output_end) = y(output_start : output_end) + real(convolution_result(8:end));
        end 
         if i==num_segments+1
        output_start =(i-1)*segment_length+1;
        output_end = output_start + length(x_keep_segment)-1;
        y(output_start : output_end) = y(output_start : output_end) + real(convolution_result(8:8+length(x_keep_segment)-1));
        end 
       
    end
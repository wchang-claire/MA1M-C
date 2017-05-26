function [rx_bit]=LDPC_HARDdecoder(rx, H, iteration)
%iteration = 10;
% Inputs: rx: Received noised bits
%         H:  Parity check matrix
%         iteration: Numbers of maximum iterations
% Output: rx_bit: Denoised bits = Received bits of the channel
%% Initializations
[M N] = size(H); % Getting the size of H
ci = rx';
rji = zeros(M, N); % Initialization of messages from check nodes
qij = H.*repmat(ci, M, 1); % Asscociate the ci matrix with non-zero elements of H, message of variable nodes
% sij = mod(ci*H',2);
%% Iterations
for n = 1:iteration
   for i = 1:M % Horizontal: check nodes -> variable nodes, each column defines a Variable
      c1 = find(H(i, :)); % Find non-zeros in the column = 1s of variable nodes for the i-th check node
      for k = 1:length(c1) % Computation of the sent-back values of variable nodes
         rji(i, c1(k)) = mod(sum(qij(i, c1)) + qij(i, c1(k)), 2);
      end % for k
%       sij = mod(rji(i,:)*H',2);
%                % If syndrome is zero
%       if sij == 0
%          n = iteration;
%       end
   end % for i


   for j = 1:N % Vertical: variable nodes -> check nodes, each row defines a Check
      f1 = find(H(:, j)); % Find non-zero in the row = 1s of check nodes (received)
      numOfOnes = length(find(rji(f1, j))); % Number of 1s in a row = number of 1s in the messages from check nodes for a given c(j) (sent)
      for k = 1:length(f1)       
         % Update qij, set '1' for majority of 1s else '0', excluding f1(k)
         if numOfOnes + ci(j) >= length(f1) - numOfOnes + rji(f1(k), j)
            qij(f1(k), j) = 1;
         else
            qij(f1(k), j) = 0;
         end
      end % for k
          
      % Bit decoding, if ci*Ht=0 or exceeds iteration, stop
      if numOfOnes + ci(j) >= length(f1) - numOfOnes
         rx_bit(j) = 1;
      else
         rx_bit(j) = 0;
      end
             
   end % for j
end % for n
%end
% %% Extraction of data bits
% if rx_bit(1:M)==rx_bit(M+1:2*M)
%     index=1;
% else
%     index=0;
% end
% decoded_bit = zeros(bin_width, length(rx_bit(1, M+1:end)));
% for index=1:bin_width
%     decoded_bit(index,:) = rx_bit(index, M+1:end);
% end
end % for function
function out = qp_call_counter(cmd, val)
% qp_call_counter
% persistent counter utility for counting QP solves/calls
%
% usage:
%   qp_call_counter('reset');
%   qp_call_counter('inc');           % +1
%   n = qp_call_counter('get');       % returns current count
%   qp_call_counter('set', 123);      % set to 123

    persistent cnt
    if isempty(cnt), cnt = 0; end

    if nargin < 1 || isempty(cmd)
        cmd = 'get';
    end

    switch lower(cmd)
        case 'reset'
            cnt = 0;
            out = cnt;
        case 'inc'
            cnt = cnt + 1;
            out = cnt;
        case 'get'
            out = cnt;
        case 'set'
            if nargin < 2, error('qp_call_counter:set needs val'); end
            cnt = val;
            out = cnt;
        otherwise
            error('qp_call_counter: unknown cmd %s', cmd);
    end
end

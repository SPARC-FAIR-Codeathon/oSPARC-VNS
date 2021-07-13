% SEMAPHORE  Interfaces with POSIX semaphore.
%
%   This mex file provides an interface with the POSIX semaphore
%   functionality. For more information, see [1]. KEY must be a signed integer.
%   VAL must be a positive integer.
%
%   SEMAPHORE('create',KEY,VAL)
%      Initializes a semaphore which can later by accessed by KEY. The
%      argument VAL specifies the initial value for the semaphore.
%
%   SEMAPHORE('destroy',KEY)
%      Destroys the semaphore indexed by KEY. Destroying a semaphore that
%      other processes or threads are currently blocked on (in
%      'wait') produces undefined behavior. Using a semaphore that
%      has been destroyed produces undefined results, until the semaphore
%      has been reinitialized using 'init'.
%
%   SEMAPHORE('wait',KEY)
%      Decrements (locks) the semaphore indexed by KEY. If the
%      semaphore's value is greater than zero, then the decrement
%      proceeds, and the function returns, immediately. If the semaphore
%      currently has the value zero, then the call blocks until either it
%      becomes possible to perform the decrement (i.e., the semaphore
%      value rises above zero), or a signal handler interrupts the call.
%
%   SEMAPHORE('post',KEY)
%      Increments (unlocks) the semaphore indexed by KEY. If the
%      semaphore's value consequently becomes greater than zero, then
%      another process or thread blocked in a 'wait' call will be woken
%      up and proceed to lock the semaphore.
%
%   See also WHOSSHARED, SHAREDMATRIX.
%
%   Example:
%      semkey=1;
%      semaphore('create',semkey,1);
%      semaphore('wait',semkey)
%      semaphore('post',semkey)
%
%   [1] - http://en.wikipedia.org/wiki/Semaphore_(programming)
%
%   Copyright (c) 2011, Joshua V Dillon
%   Copyright (c) 2014, Andrew Smart 
%   All rights reserved.

% Joshua V. Dillon
% jvdillon (a) gmail (.) com
% Wed Aug 10 13:29:01 EDT 2011

% The semaphore documentation (from which this help is generated) is part
% of release 3.27 of the Linux man-pages project. A description of the
% project, and information about reporting bugs, can be found at
% http://www.kernel.org/doc/man-pages/.

% Copyright (c) 2013, Andrew Smart
% Copyright (c) 2011, Joshua Dillon
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Georgia Institute of Technology nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


        function x = medfiltAO(x)
            %make monotonic
            [maxval, maxind] = max(x);
            [minval, minind] = min(x);
            L = length(x);
            i = minind +1;
            v = minind;
            while i~=maxind
                ni = mod(i,L)+1;
                if x(ni)<x(i)
                    x(ni) = x(i);
                end
                i = ni;
            end
            
            i = maxind+1;
            while i~=minind
                ni = mod(i,L)+1;
                if x(ni)>x(i)
                    x(ni) = x(i);
                end
                i = ni;
            end
        end
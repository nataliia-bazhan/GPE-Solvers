function [l, E] = getDimless(obj)
l = sqrt(obj.config.hbar/(obj.config.M*obj.wr));
E = obj.config.hbar*obj.wr;
end


function [ fz, fzbar ] = decompose_df( F, V, f )
% Decompose df into holomorphic and antiholormorphic
% where fz = (fx - i fy) / 2.
df = compute_df(F, V, f);

fz = holomorphic_component(df);
fzbar = antiholomorphic_component(df);
end

function [fz] = holomorphic_component(df)
fz = (df(:,1) - 1i * df(:,2)) / 2;
end

function [fz] = antiholomorphic_component(df)
fz = (df(:,1) + 1i * df(:,2)) / 2;
end
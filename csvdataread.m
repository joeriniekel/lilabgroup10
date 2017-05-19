function [ a ] = csvdataread( max_t, row )
%CSVDATAREAD Summary of this function goes here
%   Detailed explanation goes here
% max_t = 498; % 499-1
max_t = max_t - 1;
% row = 3; %breathing f, conditie 1
a = csvread('data/bb_hr.csv',row-1,0,[0,0,row-1,max_t]);
a(:,1) = [];
a = a(row,:);
end


classdef CHL
    %CHL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        pos_tar
        pos_ant
        pos_ref
        pos_img
        chl_info 
        str_vec = @(r) exp(1j*2*pi*r);
    end
    
    methods
        function obj = CHL(pos_tar_input, pos_ant_input, pos_ref_input,K_a_half,c,fc,Num_ref)
            %CHL 构造此类的实例
            %   此处显示详细说明
            obj.pos_tar = pos_tar_input;
            obj.pos_ant = pos_ant_input;
            obj.pos_ref = pos_ref_input;
            obj.pos_img = [1, 1 ,-1].*obj.pos_tar;
            obj.chl_info = obj.chl_gen(K_a_half,c,fc,Num_ref);
        end
        
        function chl_info = chl_gen(obj,K_a_half,c,fc,Num_ref)
            % cell：阵列响应矩阵，水平角，俯仰角，时延
            diff_rec_tar = obj.pos_tar - obj.pos_ant;
            azi_tar = atan(diff_rec_tar(:,2)./diff_rec_tar(:,1));
            ele_tar = atan(diff_rec_tar(:,3)./sqrt((diff_rec_tar(:,2).^2+diff_rec_tar(:,1).^2)));
            diff_dis_tar = sqrt(sum(diff_rec_tar.^2,2))/c*fc;
            Ant_rec_tar = obj.str_vec(diff_dis_tar-min(diff_dis_tar(K_a_half+1)));
            
            diff_rec_img = obj.pos_img - obj.pos_ant;
            azi_img = atan(diff_rec_img(:,2)./diff_rec_img(:,1));
            ele_img = atan(diff_rec_img(:,3)./sqrt((diff_rec_img(:,2).^2+diff_rec_img(:,1).^2)));
            diff_dis_img = sqrt(sum(diff_rec_img.^2,2))/c*fc;
            Ant_rec_img = obj.str_vec(diff_dis_img-min(diff_dis_img(K_a_half+1)));

            diff_tar_rec_ref = zeros([numel(obj.pos_ant(:,1)),3,Num_ref]);
            azi_ref = zeros([numel(obj.pos_ant(:,1)),Num_ref]);
            ele_ref = zeros([numel(obj.pos_ant(:,1)),Num_ref]);
            diff_dis_ref = zeros([numel(obj.pos_ant(:,1)),Num_ref]);
            Ant_rec_ref = zeros([numel(obj.pos_ant(:,1)),Num_ref]);
            
            for i = 1 : Num_ref
                diff_tar_rec_ref(:,:,i) = obj.pos_ref(1,:,i) - obj.pos_ant;
                diff_dis_tar_ref = sqrt(sum((obj.pos_ref(1,:,i) - obj.pos_tar).^2,2));
                azi_ref(:,i) = atan(diff_tar_rec_ref(:,2,i)./diff_tar_rec_ref(:,1,i));
                ele_ref(:,i) = atan(diff_tar_rec_ref(:,3,i)./sqrt((diff_tar_rec_ref(:,2,i).^2+diff_tar_rec_ref(:,1,i).^2)));
                diff_dis_ref(:,i) = (sqrt(sum(diff_tar_rec_ref(:,:,i).^2,2))+diff_dis_tar_ref )/c*fc;
                Ant_rec_ref(:,i) = obj.str_vec(diff_dis_ref(:,i)-min(diff_dis_ref(K_a_half+1)));
            end

            chl_info = {Ant_rec_tar, azi_tar, ele_tar;...
                Ant_rec_img, azi_img, ele_img;...
                Ant_rec_ref, azi_ref, ele_ref;...
                };
        end
    end
end


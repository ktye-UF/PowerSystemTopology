function ctime_test = ctime_eval_model(model_type, model, x_test)
switch model_type
    case 'ppgp'
        tic
        % for i=1:10
        y_pred = predict_ppgasp_isotropic(model,x_test);
        % end
        ctime_test = toc;
        % ctime_test = ctime_test / 10
    case 'msgp'
        tic
        % for i=1:10
        y_pred_gp = uq_evalModel(myKrig, x_test);
        % end
        ctime_test = toc;
        % ctime_test = ctime_test / 10
    case 'sgp'
        tic
        % for i=1:10
%         pyrunfile("xsin.py");
        pyrunfile("sgp_predict_ctime_eval.py");
        % end
        ctime_test = toc;
        % ctime_test = ctime_test / 10
    otherwise
        error('model type not found')
end






































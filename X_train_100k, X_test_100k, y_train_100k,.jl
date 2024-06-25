X_train_100k, X_test_100k, y_train_100k, y_test_100k = train_test_split(lis_train_100k, logit_(train_target_un_100k), test_size=0.2)
model=model_it(X_train_100k, X_test_100k, y_train_100k, y_test_100k,lis_train_100k, logit_(train_target_un_100k),1,'100k')


X_train_100k, X_test_100k, y_train_100k, y_test_100k = train_test_split(lis_train_100k, logit_(nTa), test_size=0.2)
model1=model_it(X_train_100k, X_test_100k, y_train_100k, y_test_100k,lis_train_100k, logit_(nTa),0,'100k')